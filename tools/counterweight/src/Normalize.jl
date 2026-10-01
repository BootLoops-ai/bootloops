# Normalize.jl — Lee-algorithm normalization & eps-factorization driver.
#
# Pipeline (single distinguished variable x; eps in the base field Qe):
#   step 0  Fuchsian check (RatFunc.is_fuchsian).  If not Fuchsian -> STOP (report).
#   step 1  eps^0 balancing: at each singular point (incl. infinity) the residue A_k has
#           eigenvalues a + b*eps.  In a dlog/UT system the TARGET is a == 0 for all of them
#           (UT exponents are pure multiples of eps).  Off-target integer eigenvalues a != 0 are
#           removed by balance transformations between a point carrying the +a eigenvalue and a
#           point (often infinity) carrying a compensating -a (Lee's eigenvalue normalization).
#   step 2  constant eps-decoupling: with all residues of the form  A_k = a_k0 + eps*a_k1 where
#           a_k0 is the (off-target) eps^0 part, find a CONSTANT (x-independent) gauge G(eps)
#           that removes a_k0 simultaneously, leaving A = eps*Atilde.  In the cases this engine
#           is built for (the semisimple Fuchsian dlog systems coming from maximal-cut-LS UT
#           candidates) the eps^0 part is already gauge-removable by a single constant rotation,
#           which we solve as a linear system in the unknown gauge entries.
#
# This is intentionally a *targeted* implementation of Lee's algorithm: it certifies and, where
# only a constant decoupling is needed, completes the eps-factorization.  Where a genuine
# higher-Poincare-rank or non-semisimple obstruction is present it STOPS with a precise report
# rather than guessing.

module Normalize

using Nemo
using ..RatFunc
using ..Fuchsia

export epsfactor_constant, EpsResult, try_epsfactor

struct EpsResult
    ok::Bool
    Atilde            # Qex matrix (eps-form connection / eps) or nothing
    T                 # the gauge applied (Qex matrix) or nothing
    Anew              # transformed connection A' = eps*Atilde (Qex) or nothing
    report::Dict{String,Any}
end

# --------------------------------------------------------------------------------------------
# Constant eps-factorization: find constant invertible G (over Qe) with
#     G A G^{-1} = eps * Atilde   (dG/dx = 0, so the (dT)T^{-1} term vanishes)
# This works when A is already Fuchsian with residues  A_k = a_k0 + eps a_k1  and a SINGLE
# constant similarity simultaneously strips every a_k0.  We solve for G by requiring the
# eps^0 part of  G A G^{-1}  to vanish.  Equivalent linear condition (multiply through by G):
#     [eps^0 part of A] * G  -  G * 0  ... handled by working order-by-order in eps.
#
# Practical solver: write A = A0(x) + eps*A1(x) + O(eps^2) ... but A may have eps in
# denominators of residues only through the affine eigen structure; for the intended target
# class A is a POLYNOMIAL in eps of degree 1 once written over the common dlog denominator.  We require:
#     A(eps,x) = A0(x) + eps * A1(x)    (degree <= 1 in eps, checked).
# Then  G A0 G^{-1} = 0   AND   G A1 G^{-1} = Atilde(x).
# The first condition forces A0(x) ≡ 0 after gauge.  If A0(x) is already 0 -> done (G=I).
# Otherwise we seek constant G with  G A0(x) = 0 * G  i.e.  G A0(x) G^{-1} = 0  for ALL x,
# i.e. A0(x) is conjugated to 0 — only possible if A0 ≡ 0.  So a *constant* gauge can only
# work when A0 ≡ 0 ALREADY.  The non-trivial eps^0 removal needs an x-dependent (balance)
# gauge, handled in try_epsfactor via balances; the constant routine is the final cleanup.
# --------------------------------------------------------------------------------------------

"Split A (Qex) as A0(x) + eps*A1(x) + ... ; return (degree_in_eps, [A0,A1,...] as Qex matrices)."
function eps_decompose(ctx::Ctx, A)
    eps = ctx.eps
    m, n = size(A)
    # Robust eps-Taylor decomposition.  Entries may be rational in eps (eps in denominators is
    # allowed, as long as eps=0 is not a pole — true for Fuchsian connections of eps-factorizable
    # systems).  We extract the eps^0 and eps^1 parts as Qex (eps-free) matrices via the
    # eps-series of each entry's x-power coefficients.  Returns (1, [A0, A1]); higher orders are
    # folded into the eps-form check (is_epsform) rather than enumerated here.
    A0 = zero_matrix(ctx.Qex, m, n)
    A1 = zero_matrix(ctx.Qex, m, n)
    for i in 1:m, j in 1:n
        a = A[i,j]; iszero(a) && continue
        A0[i,j] = _qex_eps_taylor(ctx, a, 0)
        A1[i,j] = _qex_eps_taylor(ctx, a, 1)
    end
    return 1, [A0, A1]
end

# eps^k Taylor coefficient of a Qex element, returned as an eps-free Qex element.
# Works when eps=0 is regular for the element (no eps pole); for Fuchsian eps-factorizable
# connections this holds.  Implementation: expand each x-power coefficient (a Qe=Q(eps) element)
# as an eps-series and read coefficient k.
function _qex_eps_taylor(ctx::Ctx, a, k::Int)
    num = numerator(a); den = denominator(a)
    Px = parent(num)
    # x-power coefficients are Qe elements; build the eps^k coefficient poly in x for num and den
    # then combine via the eps-series of num/den.  Easiest exact route: since den is the SAME for
    # all x-powers, write a = (Σ_i num_i x^i)/(Σ_i den_i x^i); but den may also carry eps.
    # We expand the whole rational function in eps using Qe-series arithmetic per x-power is hard;
    # instead, multiply out: find common eps-series.  Use Nemo's power series ring in eps.
    return _eps_series_coeff(ctx, a, k)
end

# Extract eps^k coefficient by substituting a power-series variable for eps.  We build, per
# x-power, the Qe element, take its eps-Laurent expansion, and assemble.  To keep things exact
# and simple we use the fact that A's entries are RATIONAL in eps with eps=0 regular: write
# entry = N(x,eps)/D(x,eps); the eps^k coeff is obtained by polynomial long division of the
# eps-series N * D^{-1}.  We implement this with a small dense eps-series (truncated at order
# k+ slack) over the field Qe-without-eps = Q, carrying x-coefficients as Q(x) rationals.
function _eps_series_coeff(ctx::Ctx, a, k::Int)
    # Represent a as a rational function in eps with coefficients in Q(x).  We map x-rational
    # coefficients into a univariate rational function field Qx2 = Q(x), and eps into a series.
    # Simpler: use Nemo's ability to evaluate the eps-series by clearing denominators numerically
    # in the field — we compute the k-th eps-derivative at eps=0 divided by k!.
    f = a
    for _ in 1:k
        f = _deps(ctx, f)
    end
    val = _eval_eps0(ctx, f)
    return val // ctx.Qex(factorial(big(k)))
end

# d/deps of a Qex element (eps is in the base field Qe).
function _deps(ctx::Ctx, f)
    num = numerator(f); den = denominator(f)
    Px = parent(num)
    # derivative wrt eps acts on the Qe coefficients of the x-polynomials num, den
    dnum = _poly_deps(Px, num)
    dden = _poly_deps(Px, den)
    nf = ctx.Qex(num); df = ctx.Qex(den)
    dnf = ctx.Qex(dnum); ddf = ctx.Qex(dden)
    return dnf//df - nf*ddf//(df*df)
end

# d/deps applied coefficient-wise to a poly over Qe in x
function _poly_deps(Px, p)
    Qe = base_ring(Px)
    cs = elem_type(Qe)[]
    for i in 0:degree(p)
        push!(cs, derivative(coeff(p, i)))   # derivative of a Q(eps) element wrt eps
    end
    return Px(cs)
end

# eps-valuation of a Qe = Q(eps) element (order of zero at eps=0; negative = pole order).
function _qe_epsval(c)
    cn = numerator(c); cd = denominator(c)
    vn = 0; while vn <= degree(cn) && iszero(coeff(cn, vn)); vn += 1; end
    vd = 0; while vd <= degree(cd) && iszero(coeff(cd, vd)); vd += 1; end
    return vn - vd
end

# multiply each Qe coefficient of an x-poly by s (used to clear common eps-content).
function _poly_scale_qe(Px, p, s)
    cs = elem_type(base_ring(Px))[]
    for i in 0:degree(p)
        push!(cs, coeff(p, i) * s)
    end
    return Px(cs)
end

# evaluate a Qex element at eps=0, returning an eps-free Qex element (in x only).  Requires
# eps=0 regular.  Done by evaluating each x-power Qe coefficient at eps=0.
# Subtlety: Nemo's canonical Qex representation (monic-in-x denominator)
# can push eps-poles into individual x-coefficients of num/den even when the ELEMENT is
# eps-regular (e.g. den = x + 1/eps after a 1/eps-carrying gauge).  Clear the common
# eps-content of num and den first; only a genuinely negative valuation difference is a pole.
function _eval_eps0(ctx::Ctx, f)
    iszero(f) && return f
    num = numerator(f); den = denominator(f)
    Px = parent(num)
    mn = minimum(_qe_epsval(coeff(num, i)) for i in 0:degree(num) if !iszero(coeff(num, i)))
    md = minimum(_qe_epsval(coeff(den, i)) for i in 0:degree(den) if !iszero(coeff(den, i)))
    if mn - md < 0
        error("eps=0 is a pole")
    elseif mn - md > 0
        return zero(ctx.Qex)   # regular and vanishing at eps=0
    end
    if mn != 0
        num = _poly_scale_qe(Px, num, ctx.eps^(-mn))
    end
    if md != 0
        den = _poly_scale_qe(Px, den, ctx.eps^(-md))
    end
    n0 = _poly_eval_eps0(Px, num)
    d0 = _poly_eval_eps0(Px, den)
    return ctx.Qex(n0) // ctx.Qex(d0)
end

function _poly_eval_eps0(Px, p)
    Qe = base_ring(Px)
    cs = elem_type(Qe)[]
    for i in 0:degree(p)
        c = coeff(p, i)
        cn = numerator(c); cd = denominator(c)
        d0 = coeff(cd, 0)
        iszero(d0) && error("eps=0 is a pole")
        push!(cs, Qe(coeff(cn, 0) // d0))
    end
    return Px(cs)
end

# eps-degree of a Qex element = max over x-power coeffs of (eps-numerator-degree). Requires
# the element be polynomial in eps (eps not in any denominator). Returns -1 for zero.
function _qex_eps_degree(a)
    iszero(a) && return -1
    num = numerator(a); den = denominator(a)
    dnum = _polyQe_eps_degree(num)
    dden = _polyQe_eps_degree(den)
    if dden != 0
        error("eps appears in x-denominator; not polynomial in eps: $a")
    end
    return dnum
end

# poly over Qe in x: max eps-degree of any coefficient
function _polyQe_eps_degree(p)
    md = 0
    for i in 0:degree(p)
        c = coeff(p, i)
        iszero(c) && continue
        cn = numerator(c); cd = denominator(c)
        if degree(cd) >= 1
            error("eps in denominator of coefficient")
        end
        md = max(md, degree(cn))
    end
    return md
end

# eps-Taylor coefficients of a Qex element, as Qex elements (eps-free), up to order maxd.
function _qex_eps_coeffs(ctx::Ctx, a, maxd)
    # a = num(x;eps)/den(x;eps), den eps-free.  Extract eps^k coeff of num, divide by den.
    num = numerator(a); den = denominator(a)
    Px = parent(num)
    out = Vector{Any}(undef, maxd+1)
    for k in 0:maxd
        # build poly in x whose x^i coeff is eps^k-coeff of num's x^i coeff
        ck = _polyQe_eps_coeff(Px, num, k)
        out[k+1] = ctx.Qex(ck) // ctx.Qex(den)
    end
    return out
end

# given poly p over Qe in x, return poly over Qe in x whose coeffs are the eps^k coeff (still
# Qe-valued but eps-free) of each coeff of p.
function _polyQe_eps_coeff(Px, p, k)
    Qe = base_ring(Px)
    coeffs = elem_type(Qe)[]
    for i in 0:degree(p)
        c = coeff(p, i)        # Qe element (Q(eps))
        cn = numerator(c); cd = denominator(c)
        # cd eps-free (checked); eps^k coeff of cn / cd
        v = (degree(cn) >= k ? coeff(cn, k) : QQ(0)) // coeff(cd, 0)
        push!(coeffs, Qe(v))
    end
    return Px(coeffs)
end

"""
    try_epsfactor(ctx, A; balances_pts=nothing, max_balance_rounds=200) -> EpsResult

Top-level: certify, and reduce to eps-form via (a) eigenvalue balancing to clear eps^0
eigenvalues, then (b) verify A0 ≡ 0 so A = eps*Atilde.  If already eps-form, returns it.
Records a full report.  STOPS (ok=false) with diagnostics if it cannot reach eps-form by the
implemented moves.
"""
function try_epsfactor(ctx::Ctx, A; max_balance_rounds::Int=400, verbose::Bool=false,
                       moser_max_iter::Int=200)
    report = Dict{String,Any}()
    report["fuchsian_input"] = is_fuchsian(ctx, A)
    Tpre = identity_matrix(ctx.Qex, size(A,1))
    if !report["fuchsian_input"]
        # Stage (1): Barkatou–Moser pole-order reduction.  Previously this was a hard STOP;
        # IBP/Kira-derived connections (e.g. C₃) routinely arrive with high-order apparent
        # poles that are removable by a rational gauge.
        mr = Fuchsia.moser_reduce(ctx, A; max_iter=moser_max_iter, verbose=verbose)
        report["moser"] = Dict("ok"=>mr.ok, "iters"=>mr.iters,
                               "profile_in"=>Fuchsia._profstr(mr.profile_in),
                               "profile_out"=>Fuchsia._profstr(mr.profile_out),
                               "stop"=>mr.stop)
        if mr.irregular !== nothing
            report["moser"]["irregular"] = mr.irregular
        end
        if !mr.ok
            report["stop"] = "input not Fuchsian and Moser reduction failed: " * mr.stop
            return EpsResult(false, nothing, mr.T, mr.A, report)
        end
        A    = mr.A
        Tpre = mr.T
        report["fuchsian_after_moser"] = is_fuchsian(ctx, A)
    end

    # Already eps-form?
    ok0, B0 = is_epsform(ctx, A)
    if ok0
        report["already_epsform"] = true
        return EpsResult(true, B0, Tpre, A, report)
    end

    Acur = A
    Ttot = Tpre
    rounds = 0
    nilgrade_rounds = 0
    decouple_rounds = 0     # cap neutral constant-decouple retries (loop guard, cf. sec169)
    history = String[]
    while rounds < max_balance_rounds
        rounds += 1
        # decompose eps^0 part
        maxd, comps = eps_decompose(ctx, Acur)
        A0 = comps[1]
        if iszero(A0)
            ok, B = is_epsform(ctx, Acur)
            if ok
                report["balance_rounds"] = rounds-1
                report["history"] = history
                report["final_certificate"] = certify_epsform(ctx, Acur)
                return EpsResult(true, B, Ttot, Acur, report)
            else
                # A0 = 0 but A/eps is still eps-dependent.  Remove the eps-dependence of A/eps with
                # a CONSTANT-in-x gauge G(eps) (fast linear-algebra path: no ODE).  If that does not
                # suffice (x-dependent regrading needed), fall back to the order-by-order Magnus
                # ODE solver.
                Gok, Greg, Areg = _constant_regrade(ctx, Acur)
                if Gok
                    push!(history, "constant-in-x eps-regrading gauge G(eps)")
                    Ttot = Greg * Ttot
                    ok2, B2 = is_epsform(ctx, Areg)
                    if ok2
                        report["balance_rounds"] = rounds-1
                        report["history"] = history
                        report["used_constant_regrade"] = true
                        report["final_certificate"] = certify_epsform(ctx, Areg)
                        return EpsResult(true, B2, Ttot, Areg, report)
                    end
                end
                Mok, Tmag, Amag = _magnus_regrade(ctx, Acur)
                if Mok
                    push!(history, "Magnus eps-regrading gauge exp(Σ eps^j W_j(x))")
                    Ttot = Tmag * Ttot
                    ok3, B3 = is_epsform(ctx, Amag)
                    if ok3
                        report["balance_rounds"] = rounds-1
                        report["history"] = history
                        report["used_magnus"] = true
                        report["final_certificate"] = certify_epsform(ctx, Amag)
                        return EpsResult(true, B3, Ttot, Amag, report)
                    else
                        Acur = Amag; continue
                    end
                end
                report["stop"] = "A0 vanished but eps-dependence of A/eps not removable by constant " *
                                 "or Magnus regrading (genuine higher-weight coupling — METHOD sec.6)"
                report["balance_rounds"] = rounds-1
                report["history"] = history
                return EpsResult(false, nothing, Ttot, Acur, report)
            end
        end
        # choose a balance to reduce ||A0|| (the eps^0 obstruction).
        moved, Acur2, T = _one_balance!(ctx, Acur, A0; verbose=verbose)
        verbose && println(stderr, "  [normalize] round $rounds: obs=",
                           _obstruction_norm(ctx, Acur), " moved=", moved,
                           moved ? " :: "*T[2] : "")
        if !moved
            # Balances cannot reduce the eps^0 obstruction (eps^0 residue eigenvalues are all
            # zero -> nilpotent / coupled eps^0 piece, OR higher-eps contamination).
            # First try Lee's "factor eps" stage: a CONSTANT-in-x gauge S(eps)⁻¹ built from the
            # Q(eps)-eigenbasis of one residue.  When all eigvals are ∝ eps this brings the
            # nilpotent eps^0 part to zero in one shot.
            Fok, Fmat, Afac = _factor_eps_eigenbasis(ctx, Acur; verbose=verbose)
            if Fok
                push!(history, "factor-eps: constant Q(eps)-eigenbasis gauge")
                Ttot = Fmat * Ttot
                ok, B = is_epsform(ctx, Afac)
                if ok
                    report["balance_rounds"] = rounds-1
                    report["history"] = history
                    report["used_factor_eps_eigenbasis"] = true
                    report["final_certificate"] = certify_epsform(ctx, Afac)
                    return EpsResult(true, B, Ttot, Afac, report)
                else
                    Acur = Afac; continue
                end
            end
            # Fallback: rational coboundary / Magnus.  Bounded retries: on ledgers where the
            # residual obstruction sits at ∞ or an algebraic orbit, the finite-pole coboundary
            # "succeeds" without global progress every round (measured on sec169) — cap it.
            Gok, Gmat, Adec = decouple_rounds <= 2*size(Acur,1) + 2 ?
                              _constant_epsdecouple(ctx, Acur) : (false, nothing, Acur)
            if Gok
                decouple_rounds += 1
                push!(history, "constant eps-decoupling gauge G(eps) (Magnus)")
                Ttot = Gmat * Ttot
                ok, B = is_epsform(ctx, Adec)
                if ok
                    report["balance_rounds"] = rounds-1
                    report["history"] = history
                    report["used_constant_gauge"] = true
                    report["final_certificate"] = certify_epsform(ctx, Adec)
                    return EpsResult(true, B, Ttot, Adec, report)
                else
                    Acur = Adec   # progress; loop again (balances may now apply)
                    continue
                end
            end
            # Lee §5 non-semisimple: nilpotent A0 with no eigenbasis / no rational coboundary.
            # Apply diagonal ε^{-L} grading (depth in the nilpotent flag of A0).  Bounded by n.
            if nilgrade_rounds < size(Acur,1)
                Nok, U, Agra, L = _nilpotent_eps_grading(ctx, A0, Acur; verbose=verbose)
                if Nok
                    nilgrade_rounds += 1
                    push!(history, "eps^{-L} nilpotent grading U=diag(eps^-L), L=$(L)")
                    Ttot = U * Ttot
                    ok, B = is_epsform(ctx, Agra)
                    if ok
                        report["balance_rounds"] = rounds-1
                        report["history"] = history
                        report["used_nilpotent_grading"] = true
                        report["nilpotent_grading_L"] = L
                        report["final_certificate"] = certify_epsform(ctx, Agra)
                        return EpsResult(true, B, Ttot, Agra, report)
                    else
                        Acur = Agra; continue   # re-enter balance loop on A_new
                    end
                end
            end
            report["stop"] = "no eps^0-reducing balance found and constant eps-decoupling failed; " *
                             "the eps^0 obstruction is not removable by the implemented moves " *
                             "(genuine coupled-DE obstruction — METHOD sec.6 / higher Poincare rank)"
            report["balance_rounds"] = rounds-1
            report["history"] = history
            report["residual_A0_nonzero"] = true
            return EpsResult(false, nothing, Ttot, Acur, report)
        end
        push!(history, T[2])
        Acur = Acur2
        Ttot = T[1] * Ttot
    end
    report["stop"] = "max balance rounds ($max_balance_rounds) exceeded"
    return EpsResult(false, nothing, Ttot, Acur, report)
end

# --------------------------------------------------------------------------------------------
# Constant-in-x eps-decoupling AND nilpotent eps^0 removal.
#
# The eps^0 part A0(x) = Σ_k N_k/(x-x_k) (N_k = eps^0 residues) obstructs eps-form.  When the
# residue *eigenvalues* are all zero (so balances do nothing) but the system is genuinely
# eps-factorizable, there is a RATIONAL gauge W(x) with  A0 = (dW/dx) W^{-1},  i.e. A0 is a
# rational coboundary; then T = W^{-1} removes A0 and the transformed system starts at eps^1.
#
# We search for W as a constant-in-eps matrix whose entries are rational in x with denominators
# drawn from the dlog letters (x - x_k).  Concretely we use the ansatz that W is a *constant*
# (x-independent) invertible matrix times a finite product of unipotent factors built from the
# nilpotent residues — equivalently we solve the linear coboundary equation order by order.
#
# Implementation (decidable & exact): seek W(x) = I + Σ_k C_k/(x-x_k) + ... is NOT polynomial in
# the right basis; instead we solve directly for the gauge by treating  dW = A0 W  as a linear
# recursion on an ansatz  W = Σ_{|α|≤L} W_α * m_α(x),  m_α monomials in the partial fractions
# 1/(x-x_k) and 1, with unknown constant matrices W_α, requiring dW - A0 W = 0 and W(∞)=I.  For
# nilpotent A0 the recursion terminates at L = nilpotency degree, giving an exact rational W.
# --------------------------------------------------------------------------------------------

# Magnus eps-regrading: A0 = 0, A = Σ_{j>=1} eps^j A_j(x).  Find gauge T(eps,x) so that
# A' = T A T^{-1} + (dT)T^{-1} = eps*A_1(x).  We solve order-by-order in eps for T = I + Σ eps^j W_j(x):
# the eps^{k} (k>=2) condition is the linear ODE
#     dW_{k-1}/dx + [W_{k-1}, A_1] = -A_k - Σ_{known lower-order brackets},
# which we solve for a RATIONAL W_{k-1}(x) using an ansatz over the dlog-letter partial fractions.
# Returns (ok, T, A').
function _magnus_regrade(ctx::Ctx, A; max_order::Int=12)
    n = size(A,1)
    # collect A_j for j = 1 .. D using the eps-Taylor extractor; find contamination degree D.
    Aj = Vector{Any}()       # Aj[k] = A_k (k>=1), as Qex matrices
    k = 1
    while true
        Ak = _eps_taylor_mat(ctx, A, k)
        push!(Aj, Ak)
        # stop when subsequent orders are zero up to a margin; detect degree by checking a few extra
        if k >= max_order
            break
        end
        # peek next two orders: if both zero, we have the full polynomial part... but A may be
        # rational in eps (infinite series).  We instead bound by max_order and verify at the end.
        k += 1
    end
    A1 = Aj[1]
    target = ctx.eps * A1
    # Build T = I + Σ_{j=1}^{m} eps^j W_j with unknown rational W_j; solve sequentially.
    Wlist = Vector{Any}()    # Wlist[j] = W_j
    poles, _ = singular_points_x(ctx, A1)
    T = identity_matrix(ctx.Qex, n)
    for j in 1:(length(Aj)-1)
        # condition at eps^{j+1}: coefficient of eps^{j+1} in (T A T^{-1} + dT T^{-1} - target) = 0.
        # Using current T (with W_1..W_{j-1} known, W_j unknown), the eps^{j+1} coefficient is
        #   dW_j/dx + [W_j, A_1] + R_j(x)  where R_j is built from known data.
        # We compute R_j by forming the residual with W_j = 0, extracting eps^{j+1}.
        Ttrial = T   # W_j not yet added
        resid = _gauge_residual(ctx, A, Ttrial, target)
        Rj = _eps_taylor_mat(ctx, resid, j+1)     # the inhomogeneous term (with W_j absent)
        # solve dW_j/dx + [W_j, A_1] = -Rj for rational W_j
        Wok, Wj = _solve_commutator_ode(ctx, A1, -Rj, poles)
        Wok || return (false, nothing, A)
        push!(Wlist, Wj)
        T = T + (ctx.eps^j) * Wj
    end
    # apply T and verify it actually reached eps-form (else report failure, do not loop)
    Aprime = apply_gauge(ctx, A, T)
    if is_epsform(ctx, Aprime)[1]
        return (true, T, Aprime)
    end
    return (false, nothing, A)
end

# eps^k Taylor coefficient of a Qex MATRIX
function _eps_taylor_mat(ctx::Ctx, A, k)
    m,n = size(A)
    B = zero_matrix(ctx.Qex, m, n)
    for i in 1:m, j in 1:n
        iszero(A[i,j]) && continue
        B[i,j] = _qex_eps_taylor(ctx, A[i,j], k)
    end
    return B
end

# residual matrix R(eps,x) = T A T^{-1} + (dT)T^{-1} - target  (as Qex matrix)
function _gauge_residual(ctx::Ctx, A, T, target)
    return apply_gauge(ctx, A, T) - target
end

# Solve  dW/dx + [W, A1] = S(x)  for a rational matrix W, using the partial-fraction ansatz.
# A1 = Σ_k a_k/(x-x_k).  We expand W = Σ_k C_k/(x-x_k) + Σ_k Σ_{m>=2} C_{k,m}/(x-x_k)^m + W_poly,
# and match.  For dlog systems with the eps-form already in A1, the solution is rational and the
# leading partial-fraction match gives  -C_k + [C_k, a_k]_residue... — we solve the linear system
# directly by collecting the coefficient of each partial-fraction basis function on both sides.
function _solve_commutator_ode(ctx::Ctx, A1, S, poles)
    n = size(A1,1)
    # Basis of rational functions for W entries: 1, and 1/(x-x_k)^m for each pole and m=1..M.
    # Determine M from the max pole order appearing in S and A1 (commutator can raise order by 1).
    M = 1 + max(_max_pole_order(ctx, S, poles), 1)
    # assemble candidate basis functions phi_b(x)
    phis = Any[ctx.Qex(1)]
    polelist = collect(poles)
    for p in polelist
        for m in 1:M
            push!(phis, ctx.Qex(1)//(ctx.x - ctx.Qex(p))^m)
        end
    end
    Nb = length(phis)
    # unknowns: W = Σ_b Wmat_b * phi_b(x), Wmat_b are n×n constant (Q) matrices => n^2 * Nb scalars.
    # Equation L(W) := dW/dx + [W,A1] - S = 0.  L is linear in the Wmat_b.  We build it as a big
    # linear system over Q by evaluating the equation at enough x-sample points (n^2 components ×
    # enough rational x values) — exact rational interpolation.
    # Number of scalar unknowns:
    nun = n*n*Nb
    # We collect equations by expanding L(W) into partial fractions and matching each basis fn &
    # matrix entry.  Simpler robust route: sample at (Nb + M + 3) distinct rational x-points and
    # require each matrix entry of L(W) to vanish there -> linear system, then solve, then VERIFY
    # the exact symbolic equation.  Build via linear algebra over QQ.
    xs = _distinct_rational_points(ctx, polelist, Nb + M + 5)
    # Precompute phi_b at each sample (as Qe rationals, eps-free => QQ via eval at any eps; phis
    # are eps-free so evaluate at eps arbitrary -> QQ).  We evaluate A1 too (eps-free since it is
    # the dlog matrix).  Build rows.
    rows = Vector{Vector{QQFieldElem}}()
    rhs = Vector{QQFieldElem}()
    # For each sample point xv and each (i,j) entry we get one equation:
    #   Σ_b [ phi_b'(xv) Wmat_b + (Wmat_b A1(xv) - A1(xv) Wmat_b)*phi_b(xv) ]_{ij} = S(xv)_{ij}
    A1v = Dict{Any,Any}()
    for xv in xs
        A1v[xv] = _eval_xv(ctx, A1, xv)
    end
    # column index for (b, r, c) unknown Wmat_b[r,c]
    colidx(b,r,c) = ((b-1)*n + (r-1))*n + c
    for xv in xs
        A1m = A1v[xv]
        Sm = _eval_xv(ctx, S, xv)
        # phi values & derivatives at xv
        phv = [_eval_scalar_xv(ctx, p, xv) for p in phis]
        dph = [_eval_scalar_xv(ctx, derivative(p), xv) for p in phis]
        for i in 1:n, j in 1:n
            row = zeros(QQ, nun)
            # contribution of each unknown Wmat_b[r,c]:
            #   d/dx term:  phi_b'(xv) * (Wmat_b)_{ij}  -> only when r=i,c=j
            #   commutator: ( Wmat_b A1 - A1 Wmat_b )_{ij} * phi_b(xv)
            #     (Wmat_b A1)_{ij} = Σ_c' Wmat_b[i,c'] A1[c',j]
            #     (A1 Wmat_b)_{ij} = Σ_r' A1[i,r'] Wmat_b[r',j]
            for b in 1:Nb
                # d/dx term
                row[colidx(b,i,j)] += dph[b]
                # commutator: + phi_b * Wmat_b[i,c'] A1[c',j]
                for cp in 1:n
                    row[colidx(b,i,cp)] += phv[b]*A1m[cp,j]
                end
                # - phi_b * A1[i,r'] Wmat_b[r',j]
                for rp in 1:n
                    row[colidx(b,rp,j)] -= phv[b]*A1m[i,rp]
                end
            end
            push!(rows, row)
            push!(rhs, Sm[i,j])
        end
    end
    # solve the (overdetermined but consistent) linear system over QQ
    Amat = matrix(QQ, length(rows), nun, reduce(vcat, [permutedims(r) for r in rows]))
    bvec = matrix(QQ, length(rhs), 1, rhs)
    sol = _solve_linear_qq(Amat, bvec)
    sol === nothing && return (false, nothing)
    # assemble W
    W = zero_matrix(ctx.Qex, n, n)
    for b in 1:Nb, r in 1:n, c in 1:n
        v = sol[colidx(b,r,c),1]
        iszero(v) && continue
        W += ctx.Qex(v) * (zero_one(ctx,n,r,c)) * phis[b]
    end
    # verify exactly: dW/dx + [W,A1] == S
    lhs = dmat_dx(ctx, W) + (W*A1 - A1*W)
    if lhs == S
        return (true, W)
    end
    return (false, nothing)
end

function zero_one(ctx::Ctx, n, r, c)
    E = zero_matrix(ctx.Qex, n, n); E[r,c] = ctx.Qex(1); return E
end

function _max_pole_order(ctx::Ctx, M, poles)
    mx = 0
    m,n = size(M)
    for i in 1:m, j in 1:n
        iszero(M[i,j]) && continue
        for p in poles
            mx = max(mx, _pole_mult(ctx, M[i,j], p))
        end
    end
    return mx
end

# distinct rational x sample points avoiding the poles
function _distinct_rational_points(ctx::Ctx, polelist, howmany)
    pvals = Set(string.(polelist))
    pts = QQFieldElem[]
    v = 2
    while length(pts) < howmany
        cand = QQ(v)
        if !(string(cand) in pvals)
            push!(pts, cand)
        end
        v += 1
        if v > howmany + length(polelist) + 50
            break
        end
    end
    return pts
end

# evaluate an eps-free Qex matrix at x = xv (QQ) -> QQ matrix
function _eval_xv(ctx::Ctx, M, xv)
    m,n = size(M)
    out = zero_matrix(QQ, m, n)
    for i in 1:m, j in 1:n
        out[i,j] = _eval_scalar_xv(ctx, M[i,j], xv)
    end
    return out
end

# evaluate an eps-free Qex scalar at x = xv (QQ) -> QQ
function _eval_scalar_xv(ctx::Ctx, f, xv)
    iszero(f) && return QQ(0)
    num = numerator(f); den = denominator(f)
    nv = _polyQe_at_x_eps0(num, xv)
    dv = _polyQe_at_x_eps0(den, xv)
    return nv // dv
end

# evaluate a poly over Qe in x at x=xv, taking eps^0 of each coeff (entries are eps-free here)
function _polyQe_at_x_eps0(p, xv)
    s = QQ(0)
    for i in 0:degree(p)
        c = coeff(p, i)
        cn = numerator(c); cd = denominator(c)
        cval = coeff(cn, 0) // coeff(cd, 0)
        s += cval * xv^i
    end
    return s
end

# solve A x = b over QQ (least-structure: use can_solve / solve)
function _solve_linear_qq(A, b)
    fl, x = can_solve_with_solution(A, b; side=:right)
    fl || return nothing
    return x
end

# Constant-in-x eps-regrading.  A0 = 0, so M := A/eps = Σ_{j>=0} eps^j M_j(x), M_0 = target dlog.
# Find a constant-in-x gauge G(eps) = I + Σ_{k>=1} eps^k G_k (G_k constant Q-matrices) so that
# G M G^{-1} is eps-free.  Because G is x-independent, the gauge action on A = eps M is pure
# conjugation: A -> eps (G M G^{-1}), Fuchsian-structure preserved, no derivative term.
# Order eps^k condition: the eps^k coefficient of G M G^{-1} must vanish (k>=1).  At order k it is
#   M_k + [G_k, M_0] + (known terms in G_1..G_{k-1}) = 0,
# a linear system for the constant matrix G_k:  G_k M_0(x) - M_0(x) G_k = -(M_k + known)(x),
# required for all x.  We solve it by sampling x at enough rational points (exact) and verify.
function _constant_regrade(ctx::Ctx, A; max_k::Int=10)
    n = size(A,1)
    eps = ctx.eps
    # M = A/eps
    M = zero_matrix(ctx.Qex, n, n)
    for i in 1:n, j in 1:n
        iszero(A[i,j]) && continue
        M[i,j] = A[i,j] // ctx.Qex(eps)
    end
    M0 = _eps_taylor_mat(ctx, M, 0)        # target dlog (eps-free)
    poles, _ = singular_points_x(ctx, M0)
    polelist = collect(poles)
    # Solve for the constant-in-x gauge G(eps) (entries in Q(eps)) directly from
    #     G M(eps,x) - M0(x) G = 0   for all x,
    # which says G conjugates M to the eps-free M0.  Linear in G over the field Q(eps); we sample
    # enough x-points and solve over Q(eps).  G is required invertible (det != 0).
    Gqe = _solve_conjugation_const(ctx, M, M0, polelist, n)
    Gqe === nothing && return (false, nothing, A)
    G = _qe_to_qex(ctx, Gqe)
    # verify and apply
    Areg = apply_gauge(ctx, A, G)
    if is_epsform(ctx, Areg)[1]
        return (true, G, Areg)
    end
    return (false, nothing, A)
end

# Solve for constant-in-x G (entries in Qe = Q(eps)) with  G*M(eps,x) = M0(x)*G  for all x.
# M, M0 are Qex matrices (M eps-dependent, M0 eps-free).  We sample x at rationals, building a
# homogeneous linear system over Qe for the n^2 unknown entries of G, then pick the solution that
# is invertible (normalize the gauge by fixing it to the identity-branch: we seek G with the same
# "shape" as a unipotent/invertible conjugator).  Returns a Qe matrix or nothing.
function _solve_conjugation_const(ctx::Ctx, M, M0, polelist, n)
    Qe = ctx.Qe
    nun = n*n
    colidx(r,c) = (r-1)*n + c
    xs = _distinct_rational_points(ctx, polelist, 3*n + 6)
    # Build homogeneous system  Σ_{r,c} coeff * G[r,c] = 0  (entries in Qe).
    rows = Vector{Vector{elem_type(Qe)}}()
    for xv in xs
        Mm  = _eval_xv_qe(ctx, M, xv)     # Qe matrix (eps-dependent)
        M0m = _eval_xv_qe(ctx, M0, xv)    # Qe matrix (eps-free values, still in Qe)
        for i in 1:n, j in 1:n
            row = [Qe(0) for _ in 1:nun]
            # (G M)_{ij} = Σ_c G[i,c] M[c,j]
            for c in 1:n
                row[colidx(i,c)] += Mm[c,j]
            end
            # -(M0 G)_{ij} = -Σ_r M0[i,r] G[r,j]
            for r in 1:n
                row[colidx(r,j)] -= M0m[i,r]
            end
            push!(rows, row)
        end
    end
    Amat = _rows_to_matrix_qe(Qe, rows, nun)
    # nullspace over Qe
    ns = _nullspace_qe(Amat)
    ns === nothing && return nothing
    (nd, basis) = ns
    nd == 0 && return nothing
    # The conjugator gauge is unique up to multiplication by matrices commuting with M0.  We pick
    # a combination from the nullspace that is invertible and reduces to a sensible G at eps=0.
    # Strategy: assemble candidate G from each nullspace vector; we want G|_{eps=0} invertible.
    # Try single basis vectors and simple combinations.
    cand = _pick_invertible_gauge(ctx, basis, n)
    cand === nothing && return nothing
    return cand
end

# evaluate Qex matrix at x=xv keeping eps (-> Qe matrix)
function _eval_xv_qe(ctx::Ctx, M, xv)
    m,n = size(M)
    Qe = ctx.Qe
    out = zero_matrix(Qe, m, n)
    for i in 1:m, j in 1:n
        iszero(M[i,j]) && continue
        out[i,j] = _eval_scalar_xv_qe(ctx, M[i,j], xv)
    end
    return out
end

function _eval_scalar_xv_qe(ctx::Ctx, f, xv)
    num = numerator(f); den = denominator(f)
    nv = _polyQe_at_x(num, xv)
    dv = _polyQe_at_x(den, xv)
    return nv // dv
end

# evaluate poly over Qe in x at x=xv (xv QQ), keeping eps -> Qe element
function _polyQe_at_x(p, xv)
    Qe = base_ring(parent(p))
    s = Qe(0)
    for i in 0:degree(p)
        s += coeff(p, i) * Qe(xv)^i
    end
    return s
end

function _rows_to_matrix_qe(Qe, rows, ncol)
    nr = length(rows)
    M = zero_matrix(Qe, nr, ncol)
    for i in 1:nr, j in 1:ncol
        M[i,j] = rows[i][j]
    end
    return M
end

# nullspace over Qe; returns (dim, basis-matrix with columns = basis vectors) or nothing
function _nullspace_qe(A)
    nd, ker = nullspace(A)
    return (nd, ker)
end

# from a nullspace basis (columns of ker, length n^2), pick a combination forming an invertible
# n×n matrix G with G|_{eps=0} invertible.  Try basis vectors and pairwise sums.
function _pick_invertible_gauge(ctx::Ctx, ker, n)
    Qe = ctx.Qe
    ncols = size(ker, 2)
    function vec_to_mat(col)
        G = zero_matrix(Qe, n, n)
        for r in 1:n, c in 1:n
            G[r,c] = ker[(r-1)*n + c, col]
        end
        return G
    end
    function ok_gauge(G)
        # invertible over Qe and finite/invertible at eps=0
        try
            isz = iszero(det(G))
            isz && return false
            G0 = eps0_part(ctx, G)
            return !iszero(det(G0))
        catch
            return false
        end
    end
    # single basis vectors
    for c in 1:ncols
        G = vec_to_mat(c)
        ok_gauge(G) && return G
    end
    # combinations Σ c_i v_i with small integer coeffs
    if ncols >= 2
        for c1 in 1:ncols, c2 in (c1+1):ncols
            for a in (1,-1,2), b in (1,-1)
                G = a*vec_to_mat(c1) + b*vec_to_mat(c2)
                ok_gauge(G) && return G
            end
        end
    end
    return nothing
end

# is a Qex matrix eps-free?
function _mat_eps_free(ctx::Ctx, M)
    m,n = size(M)
    for i in 1:m, j in 1:n
        iszero(M[i,j]) && continue
        # eps-derivative zero?
        if !iszero(_deps(ctx, M[i,j]))
            return false
        end
    end
    return true
end

# solve constant G with  G*M0(x) - M0(x)*G = R(x)  for all x, via sampling.  R, M0 eps-free Qex.
function _solve_commutator_const(ctx::Ctx, M0, R, polelist, n)
    nun = n*n
    colidx(r,c) = (r-1)*n + c
    xs = _distinct_rational_points(ctx, polelist, 2*n + 6)
    rows = Vector{Vector{QQFieldElem}}()
    rhs = QQFieldElem[]
    for xv in xs
        M0m = _eval_xv(ctx, M0, xv)
        Rm  = _eval_xv(ctx, R, xv)
        for i in 1:n, j in 1:n
            row = zeros(QQ, nun)
            # (G M0)_{ij} = Σ_c G[i,c] M0[c,j] ; (M0 G)_{ij} = Σ_r M0[i,r] G[r,j]
            for c in 1:n
                row[colidx(i,c)] += M0m[c,j]
            end
            for r in 1:n
                row[colidx(r,j)] -= M0m[i,r]
            end
            push!(rows, row)
            push!(rhs, Rm[i,j])
        end
    end
    Amat = _rows_to_matrix(rows, nun)
    bvec = matrix(QQ, length(rhs), 1, rhs)
    sol = _solve_linear_qq(Amat, bvec)
    sol === nothing && return nothing
    G = zero_matrix(QQ, n, n)
    for r in 1:n, c in 1:n
        G[r,c] = sol[colidx(r,c),1]
    end
    # verify symbolic commutator
    Gx = _qe_to_qex(ctx, G)
    if (Gx*M0 - M0*Gx) == R
        return G
    end
    return nothing
end

function _rows_to_matrix(rows, ncol)
    nr = length(rows)
    M = zero_matrix(QQ, nr, ncol)
    for i in 1:nr, j in 1:ncol
        M[i,j] = rows[i][j]
    end
    return M
end

# Qe (QQ) matrix -> Qex matrix
function _qe_to_qex(ctx::Ctx, G)
    m,n = size(G)
    out = zero_matrix(ctx.Qex, m, n)
    for i in 1:m, j in 1:n
        out[i,j] = ctx.Qex(G[i,j])
    end
    return out
end

"""
    _factor_eps_eigenbasis(ctx, A) -> (ok, G, A')

Lee "factor eps" (§5): when every residue eigenvalue is a pure multiple of eps (eps^0 part = 0
at every pole), a CONSTANT-in-x gauge G(eps) = S(eps)⁻¹ with S the Q(eps)-eigenbasis of one
residue R(p₀) brings A to eps-form.  We try each pole in turn (the one with the most distinct
eigvals first) and return the first G for which is_epsform(G A G⁻¹) holds.  G is generically
singular at eps=0 (det S = O(eps^k)); that is the expected weight-grading of the canonical UT
basis and is harmless — only the connection A' = G A G⁻¹ needs to be eps-regular.
"""
function _factor_eps_eigenbasis(ctx::Ctx, A; verbose::Bool=false)
    n = size(A,1)
    pts, _ = singular_points_x(ctx, A)
    allpts = vcat(collect(pts), [:inf])
    # quick gate: every eigval at every pole must have eps^0 part 0
    cand = Any[]
    for p in allpts
        R = (p === :inf) ? residue_at_infinity(ctx, A) : residue_matrix(ctx, A, p)
        evs = _qe_linear_eigvals(ctx, R)
        for (_, a, _) in evs
            a == 0 || return (false, nothing, A)
        end
        push!(cand, (p, R, length(evs)))
    end
    # try poles with the most DISTINCT eigvals first (best chance of a full eigenbasis)
    sort!(cand; by = t -> -t[3])
    for (p, R, _) in cand
        S = zero_matrix(ctx.Qe, n, n)
        col = 0
        for (lam, _, _) in _qe_linear_eigvals(ctx, R)
            for u in _qe_right_eigvecs(ctx, R, lam)
                col >= n && break
                col += 1
                for i in 1:n; S[i,col] = u[i,1]; end
            end
        end
        col == n || continue
        iszero(det(S)) && continue
        G = Fuchsia._qe_to_qex(ctx, inv(S))
        Anew = apply_gauge(ctx, A, G)
        if is_epsform(ctx, Anew)[1]
            verbose && println(stderr, "  [normalize] factor-eps via Q(eps)-eigenbasis at p=", p)
            return (true, G, Anew)
        end
    end
    return (false, nothing, A)
end

"""
    _nilpotent_eps_grading(ctx, A0, A) -> (ok, U, A', L)

Lee §5 non-semisimple "factor ε": when every ε^0 residue is NILPOTENT (so balances do
nothing and `_factor_eps_eigenbasis` finds no full eigenbasis), apply the diagonal weight
gauge U = diag(ε^{-L_i}), L_i = depth of row i in the nilpotent flag of A0.  Then
(U A U⁻¹)[i,j] = ε^{L_j-L_i}·A[i,j], pushing every nonzero A0[i,j] (which has L_j > L_i by
construction) to ε^{≥1}.  U is x-constant, so dU·U⁻¹ = 0.  Returns ok=false if A0's support
graph has a cycle or a nonzero diagonal (no diagonal grading exists).
"""
function _nilpotent_eps_grading(ctx::Ctx, A0, A; verbose::Bool=false)
    n = size(A,1)
    iszero(A0) && return (false, nothing, A, Int[])
    # gate: every ε^0 residue (incl. ∞) must be nilpotent
    pts, _ = singular_points_x(ctx, A0)
    for p in vcat(collect(pts), Any[:inf])
        N = (p === :inf) ? eps0_part(ctx, residue_at_infinity(ctx, A0)) :
                           eps0_part(ctx, residue_matrix(ctx, A0, p))
        _is_nilpotent(N) || return (false, nothing, A, Int[])
    end
    # incidence DAG i→j  ⇔  A0[i,j] ≠ 0  (any x); diagonal must vanish
    adj = [Int[] for _ in 1:n]; indeg = zeros(Int, n)
    for i in 1:n, j in 1:n
        iszero(A0[i,j]) && continue
        i == j && return (false, nothing, A, Int[])
        push!(adj[i], j); indeg[j] += 1
    end
    # Kahn longest-path layering: L_j = 1 + max_{i: i→j} L_i
    L = zeros(Int, n); deg = copy(indeg)
    queue = [i for i in 1:n if deg[i] == 0]; seen = 0
    while !isempty(queue)
        i = popfirst!(queue); seen += 1
        for j in adj[i]
            L[j] = max(L[j], L[i]+1)
            (deg[j] -= 1) == 0 && push!(queue, j)
        end
    end
    (seen == n && maximum(L) > 0) || return (false, nothing, A, Int[])  # cycle / no edges
    # U = diag(ε^{-L_i}); apply as similarity (x-const ⇒ no dU term)
    eps = ctx.Qex(ctx.eps)
    U = zero_matrix(ctx.Qex, n, n); Ui = zero_matrix(ctx.Qex, n, n)
    for i in 1:n; U[i,i] = inv(eps^L[i]); Ui[i,i] = eps^L[i]; end
    Anew = U * A * Ui
    verbose && println(stderr, "  [normalize] eps^{-L} nilpotent grading L=", L)
    return (true, U, Anew, L)
end

# ==============================================================================================
# Case A'' — conjugate-pair (Galois-orbit) Lee balances at irreducible algebraic letters.
#
# The rational-point Lee moves above are blind to residues sitting at the roots of an
# IRREDUCIBLE (ε-free) denominator factor q(x) of degree d ≥ 2 (e.g. a quadratic
# 6x²+17x+46 with λ = −1−2ε doubled at each conjugate root, so Σa over the rational poles alone
# is ≠ 0 and no sequence of rational balances can reach ε-form).  Same tool-blindness class
# as the R∞ case (quadratic letters at infinity vs quadratic/cubic poles at finite x).
#
# Cure: balance the WHOLE Galois orbit at once with the ℚ(ε)-rational gauge
#     T = (I − P) + f(x)·P ,     f = q̂(x) / Π_{b ∈ D} (x − b)          (raise a by +1 at
#                                                                        EVERY root of q̂)
#     (or its mirror f = Π (x−b)/q̂ to lower at the roots),
# where q̂ is the monic form of q and D is a "dump" multiset of |D| = d rational points
# (∞ absorbing the remaining degree).  P must be the ℚ(ε)-rational rank-1 projector
#     P = u vᵀ / (vᵀ u)
# built WITHOUT ever leaving ℚ(ε):
#   * u spans the JOINT rational eigenspace of the orbit-residue R(y) = Σ_j y^j R_j
#     (y ≡ a root of q̂, R_j ∈ ℚ(ε)-matrices):  R₀ u = λ u,  R_j u = 0 (j ≥ 1),
#     so R(r)u = λu at EVERY conjugate root simultaneously.  On this rational structure
#     P coincides with the spectral-pair combination P_r + P_r̄ (its ℚ-rational trace part),
#     is exactly idempotent and rank 1.
#   * vᵀ is a simultaneous LEFT eigenvector of the full residue at every dump point of D
#     (Lee's invariance condition at the lowered points), incl. R∞ when ∞ ∈ D; when ∞
#     carries multiplicity ≥ 2 in D the second-order condition vᵀ·A₂·(I−P) = 0 (A₂ = 1/x²
#     Laurent coefficient at ∞) is additionally enforced by an eigen-refinement inside the
#     R∞ eigenspace.
# Residue bookkeeping of the move: a → a+1 on the P-line at each of the d roots and
# a → a−1 at each dump point (counted with multiplicity, −(d−|D_fin|) at ∞); Fuchsianity
# conditions at the roots ((I−P)·R(root)·P = 0) hold automatically by construction of u.
# The mirror move swaps the roles of u and v.  Every candidate is hard-gated on
# Fuchsianity + ε-regularity + strict drop of the algebraic-aware obstruction norm.
# ==============================================================================================

const _ORBIT_MAXDEG = 3          # support quadratic + cubic letters
const _ORBIT_MAXCAND = 4000      # cap on (config × eigcombo) scans per round

# Collect irreducible ε-free denominator factors of degree 2.._ORBIT_MAXDEG over all entries,
# returned as monic QQ coefficient vectors cs = [c0, c1, ..., c_{d-1}, 1].
function _alg_orbit_factors(ctx::Ctx, A)
    seen = Dict{String,Vector{QQFieldElem}}()
    m, n = size(A)
    for i in 1:m, j in 1:n
        a = A[i,j]; iszero(a) && continue
        den = denominator(a)
        is_constant(den) && continue
        for (p, _) in factor(den)
            dp = degree(p)
            (2 <= dp <= _ORBIT_MAXDEG) || continue
            lc = coeff(p, dp)
            cs = QQFieldElem[]
            ok = true
            for k in 0:dp
                c = coeff(p, k) // lc                    # Qe element
                cn = numerator(c); cd = denominator(c)
                if degree(cn) >= 1 || degree(cd) >= 1
                    ok = false; break                    # ε-dependent letter: skip
                end
                push!(cs, coeff(cn, 0) // coeff(cd, 0))
            end
            ok || continue
            key = join(string.(cs), ",")
            haskey(seen, key) || (seen[key] = cs)
        end
    end
    return collect(values(seen))
end

# monic q̂(x) as a Qe[x] polynomial from its coefficient vector
function _orbit_poly(ctx::Ctx, cs)
    Px = parent(numerator(ctx.x))
    xx = gen(Px)
    return sum(Px(ctx.Qe(cs[k+1])) * xx^k for k in 0:(length(cs)-1))
end

# Residue components on the orbit of monic irreducible q̂:  R(y) = Σ_{j=0}^{d-1} y^j R_j with
# y a root of q̂.  Exact, entirely over ℚ(ε):  per entry num/den with den = q̂·d̃ (mult 1),
# residue(y) = num(y)·d̃(y)⁻¹·q̂'(y)⁻¹ computed via gcdx-inverses mod q̂.
# Returns Vector{Qe-matrix} of length d, or nothing (multiplicity ≥ 2 on the orbit).
function _orbit_residue_components(ctx::Ctx, A, cs)
    d = length(cs) - 1
    n = size(A,1)
    qpoly = _orbit_poly(ctx, cs)
    dq = derivative(qpoly)
    g, invdq, _ = gcdx(dq, qpoly)                 # dq·invdq ≡ g mod q̂, g = gcd (monic)
    isone(g) || return nothing                    # q̂ not squarefree (cannot happen: irreducible)
    Rj = [zero_matrix(ctx.Qe, n, n) for _ in 1:d]
    any_hit = false
    for i in 1:n, jj in 1:n
        f = A[i,jj]; iszero(f) && continue
        num = numerator(f); den = denominator(f)
        quo, rem1 = divrem(den, qpoly)
        iszero(rem1) || continue                  # this entry has no pole on the orbit
        rem2 = mod(quo, qpoly)
        iszero(rem2) && return nothing            # multiplicity ≥ 2: not Fuchsian on the orbit
        g2, invquo, _ = gcdx(rem2, qpoly)
        isone(g2) || return nothing
        r = mod(mod(num, qpoly) * invquo * invdq, qpoly)
        for k in 0:degree(r)
            Rj[k+1][i,jj] += coeff(r, k)
        end
        any_hit = true
    end
    any_hit || return nothing
    return Rj
end

# (d·n)×(d·n) ℚ(ε)-matrix of "multiply by R(y)" on K^n, K = ℚ(ε)[y]/(q̂): used only for the
# algebraic-aware obstruction norm (its rational eigenvalues are those of R at the roots,
# counted over the whole orbit).
function _orbit_hat_matrix(ctx::Ctx, Rj, cs)
    d = length(Rj); n = size(Rj[1], 1)
    Qe = ctx.Qe
    # reduction table: y^m (m = 0..2d-2) in the basis y^0..y^{d-1}
    ytab = Vector{Vector{QQFieldElem}}()
    for m in 0:(2d-2)
        if m < d
            v = [QQ(k == m ? 1 : 0) for k in 0:(d-1)]
        else
            prev = ytab[m]                        # y^{m-1}
            v = zeros(QQ, d)
            for k in 0:(d-2); v[k+2] = prev[k+1]; end
            top = prev[d]
            if !iszero(top)
                for k in 0:(d-1); v[k+1] -= top * cs[k+1]; end
            end
        end
        push!(ytab, v)
    end
    H = zero_matrix(Qe, d*n, d*n)
    for kin in 0:(d-1), j in 0:(d-1)
        red = ytab[j + kin + 1]
        Rjm = Rj[j+1]
        for l in 0:(d-1)
            c = red[l+1]; iszero(c) && continue
            for r in 1:n, cc in 1:n
                iszero(Rjm[r,cc]) && continue
                H[l*n + r, kin*n + cc] += Qe(c) * Rjm[r,cc]
            end
        end
    end
    return H
end

# Joint rational right eigenspace of the orbit residue: {u : R₀u = λu, R_j u = 0 (j≥1)}.
# Returns the list of Qe column vectors (ε-cleared).  Pass transposed R_j for the left space.
function _orbit_joint_eigvecs(ctx::Ctx, Rj, lam)
    n = size(Rj[1], 1)
    M = Rj[1] - lam * identity_matrix(ctx.Qe, n)
    for j in 2:length(Rj)
        M = vcat(M, Rj[j])
    end
    k, K = nullspace(M)
    return [_qe_clear_denoms(ctx, K[:, j:j]) for j in 1:k]
end

# x⁻² Laurent coefficient of A at infinity (A must be O(1/x) there).
function _laurent2_at_infinity(ctx::Ctx, A)
    m, n = size(A)
    Qe = ctx.Qe
    B = zero_matrix(Qe, m, n)
    for i in 1:m, j in 1:n
        f = A[i,j]; iszero(f) && continue
        g = subst_x_inv(ctx, f)                    # g(w) = f(1/w) = f₁w + f₂w² + ...
        h = g // ctx.x                             # f₁ + f₂w + ...
        B[i,j] = RatFunc._eval_qex_at(ctx, derivative(h), Qe(0))
    end
    return B
end

# Left eigen-refinement inside a row space: given rows W (k×n over Qe) and M (n×n), find
# combinations vᵀ ∈ rowspace(W) with vᵀM = ν vᵀ (ν ∈ Qe).  Handles k = 1, 2 exactly;
# for k > 2 falls back to testing single rows.  Returns list of 1×n row vectors.
function _refine_left_eig_in_rowspace(ctx::Ctx, W, M)
    Qe = ctx.Qe
    k = nrows(W); n = ncols(W)
    out = Any[]
    WM = W * M
    if k == 1
        # need WM ∥ W: check all 2×2 minors of [W; WM]
        ok = true
        for c1 in 1:n, c2 in (c1+1):n
            iszero(W[1,c1]*WM[1,c2] - W[1,c2]*WM[1,c1]) || (ok = false; break)
        end
        ok && push!(out, W)
        return out
    end
    if k == 2
        # pencil: find ν with rank(WM − νW) ≤ 1 via 2×2 column minors (polys in ν, deg ≤ 2)
        Rt, t = polynomial_ring(Qe, "t")
        gpol = zero(Rt)
        for c1 in 1:n, c2 in (c1+1):n
            m11 = WM[1,c1] - t*W[1,c1]; m12 = WM[1,c2] - t*W[1,c2]
            m21 = WM[2,c1] - t*W[2,c1]; m22 = WM[2,c2] - t*W[2,c2]
            dpol = m11*m22 - m12*m21
            gpol = iszero(gpol) ? dpol : gcd(gpol, dpol)
            isone(gpol) && return out              # no common ν
        end
        iszero(gpol) && (gpol = zero(Rt))          # rank ≤ 1 for ALL ν (degenerate): ν free
        nus = Any[]
        if iszero(gpol)
            push!(nus, Qe(0))
        else
            for (fac, _) in factor(gpol)
                degree(fac) == 1 || continue
                push!(nus, -coeff(fac,0)//coeff(fac,1))
            end
        end
        for nu in nus
            D = WM - nu*W
            kk, KK = nullspace(transpose(D))       # c with cᵀD = 0 ⇔ Dᵀc = 0
            for jj in 1:kk
                cvec = KK[:, jj:jj]
                v = transpose(cvec) * W            # 1×n
                iszero(v) || push!(out, _qe_clear_denoms(ctx, v))
            end
        end
        return out
    end
    # k > 2: single-row fallback
    for r in 1:k
        Wr = W[r:r, :]
        append!(out, _refine_left_eig_in_rowspace(ctx, Wr, M))
    end
    return out
end

# Algebraic-aware obstruction norm: the rational-point norm PLUS Σ|a| of the rational
# ε⁰-eigenvalues of the orbit residues over every irreducible ε-free algebraic letter
# (counted over the whole orbit via the hat model).  Used ONLY to gate Case A'' moves,
# so the behaviour of the pre-existing move ladder is untouched.
function _obstruction_norm_alg(ctx::Ctx, A)
    tot = _obstruction_norm(ctx, A)
    for cs in _alg_orbit_factors(ctx, A)
        Rj = _orbit_residue_components(ctx, A, cs)
        Rj === nothing && continue
        H0 = eps0_part(ctx, _orbit_hat_matrix(ctx, Rj, cs))
        for r in rational_roots_with_mult(charpoly_q(H0))
            tot += abs(r)
        end
    end
    return tot
end

# The orbit balance gauge: T = (I−P) + f·P with f = q̂ / Π_{b∈D_fin}(x−b)  (raise at roots)
# or f = Π_{b∈D_fin}(x−b) / q̂  (lower at roots).
function _orbit_gauge(ctx::Ctx, P, cs, dumps; lower_at_roots::Bool=false)
    n = size(P,1)
    x = ctx.x
    qx = sum(ctx.Qex(cs[k+1]) * x^k for k in 0:(length(cs)-1))
    dpoly = ctx.Qex(1)
    for b in dumps
        b === :inf && continue
        dpoly *= (x - ctx.Qex(b))
    end
    f = lower_at_roots ? (dpoly // qx) : (qx // dpoly)
    Pq = Fuchsia._qe_to_qex(ctx, P)
    Iq = identity_matrix(ctx.Qex, n)
    return (Iq - Pq) + f * Pq
end

# Enumerate dump configurations D (multisets of size d over the rational poles ∪ {∞};
# finite points distinct, ∞ with multiplicity ≤ 2 — the second-order ∞ condition is
# handled, third order is not).
function _orbit_dump_configs(polelist, d)
    configs = Vector{Any}[]
    fins = collect(polelist)
    nf = length(fins)
    function rec(start, chosen)
        nin = d - length(chosen)                   # remaining slots -> ∞ multiplicity
        if nin <= 2
            push!(configs, vcat(chosen, fill(:inf, nin)))
        end
        length(chosen) == d && return
        for idx in start:nf
            rec(idx+1, vcat(chosen, Any[fins[idx]]))
        end
    end
    rec(1, Any[])
    return configs
end

# --- K = ℚ(ε)[y]/(q̂) scalar kernel for d = 2 (pairs (a₀,a₁) ≡ a₀ + y·a₁, y² = −p·y − s) ---
_k2_add(a, b) = (a[1] + b[1], a[2] + b[2])
_k2_sub(a, b) = (a[1] - b[1], a[2] - b[2])
_k2_mul(p, s, a, b) = (a[1]*b[1] - s*(a[2]*b[2]), a[1]*b[2] + a[2]*b[1] - p*(a[2]*b[2]))
# inverse via the field norm N(a) = a₀² − p·a₀a₁ + s·a₁² (a·conj(a) = N, conj(a) = a₀−p·a₁ − y·a₁)
function _k2_inv(p, s, a)
    Nrm = a[1]^2 - p*a[1]*a[2] + s*a[2]^2
    iszero(Nrm) && return nothing
    return ((a[1] - p*a[2]) * inv(Nrm), -a[2] * inv(Nrm))
end
# reduce an x-polynomial over Qe mod q̂ to its K components (evaluation "at the root y")
_k2_ofpoly(qpoly, P) = (r = mod(P, qpoly); (coeff(r, 0), coeff(r, 1)))

# Constant ((x−root)⁰) Laurent coefficient of A on the orbit of monic irreducible q̂ (d = 2),
# in K components: B(y) = B₀ + y·B₁ with B(r) = lim_{x→r} (A − R(r)/(x−r)) at each root r.
# Exact, entirely over ℚ(ε), via mod-q̂ arithmetic.  For an entry f = num/(q̂·d̃) (orbit
# multiplicity 1): f = g(x)/(x−r) with g = num/(d̃·(x−r̄)), so the constant term is
#     g'(r) = [num'·d̃·q̂' − num·(d̃'·q̂' + d̃)] / (d̃·q̂')²   evaluated in K   (q̂'(r) = r − r̄);
# an entry with no orbit pole contributes its plain K-evaluation num(y)/den(y).  Returns
# (B₀, B₁) over Qe, or nothing (orbit multiplicity ≥ 2: not Fuchsian on the orbit).
function _orbit_const_term(ctx::Ctx, A, cs)
    n = size(A, 1)
    Qe = ctx.Qe
    p = Qe(cs[2]); s = Qe(cs[1])
    qpoly = _orbit_poly(ctx, cs)
    qp = (p, Qe(2))                                # q̂'(y) = 2y + p
    B0 = zero_matrix(Qe, n, n); B1 = zero_matrix(Qe, n, n)
    for i in 1:n, j in 1:n
        f = A[i,j]; iszero(f) && continue
        num = numerator(f); den = denominator(f)
        quo, rem1 = divrem(den, qpoly)
        local val
        if !iszero(rem1)
            dv = _k2_inv(p, s, _k2_ofpoly(qpoly, den))
            dv === nothing && return nothing       # q̂ | den after all (cannot happen)
            val = _k2_mul(p, s, _k2_ofpoly(qpoly, num), dv)
        else
            iszero(mod(quo, qpoly)) && return nothing   # multiplicity ≥ 2 on the orbit
            nm  = _k2_ofpoly(qpoly, num)
            nmd = _k2_ofpoly(qpoly, derivative(num))
            dt  = _k2_ofpoly(qpoly, quo)
            dtd = _k2_ofpoly(qpoly, derivative(quo))
            e   = _k2_mul(p, s, dt, qp)                 # d̃(y)·q̂'(y) ≠ 0 in the field K
            einv2 = _k2_inv(p, s, _k2_mul(p, s, e, e))
            einv2 === nothing && return nothing
            t1 = _k2_mul(p, s, nmd, e)
            t2 = _k2_mul(p, s, nm, _k2_add(_k2_mul(p, s, dtd, qp), dt))
            val = _k2_mul(p, s, _k2_sub(t1, t2), einv2)
        end
        B0[i,j] = val[1]; B1[i,j] = val[2]
    end
    return (B0, B1)
end

# Tr-pair projector (task-spec construction, d = 2 only): P = P_r + P_r̄ with rank-1
# spectral picks at the two conjugate roots, cross-orthogonalized so that P is idempotent:
#   u from null(R̂ − λ) (a K-eigvector u(y) = u₀ + y·u₁ of the orbit residue),
#   v from null(R̂ᵀ − λ) subject to  v(y)ᵀ·u(ȳ) = 0  in K  (kills P_r·P_r̄ cross terms),
#   P = Tr_{K/ℚ(ε)}( u(y)v(y)ᵀ / (v(y)ᵀu(y)) ) = 2Π₀ − p·Π₁  — ℚ(ε)-rational, rank 2.
# Returns a list of candidate P matrices (Qe).  cs = [s, p, 1] (monic x² + p·x + s).
function _orbit_trpair_projectors(ctx::Ctx, Rj, cs, lam)
    n = size(Rj[1], 1)
    Qe = ctx.Qe
    p = Qe(cs[2]); s = Qe(cs[1])
    Rhat  = _orbit_hat_matrix(ctx, Rj, cs)
    RhatT = _orbit_hat_matrix(ctx, [transpose(R) for R in Rj], cs)
    ku, KU = nullspace(Rhat  - lam*identity_matrix(Qe, 2n))
    kv, KV = nullspace(RhatT - lam*identity_matrix(Qe, 2n))
    (ku == 0 || kv == 0) && return Any[]
    out = Any[]
    for cu in 1:ku
        U = _qe_clear_denoms(ctx, KU[:, cu:cu])
        u0 = U[1:n, 1:1]; u1 = U[(n+1):2n, 1:1]
        # conjugate: ū = (u0 − p·u1) − y·u1
        ub0 = u0 - p*u1; ub1 = -u1
        # cross condition on v = Σ_j c_j KV_j :  v(y)ᵀ·ū(y) = 0 in K (two Qe-linear equations)
        C = zero_matrix(Qe, 2, kv)
        for j in 1:kv
            w0 = KV[1:n, j:j]; w1 = KV[(n+1):2n, j:j]
            # (w0 + y w1)ᵀ(ub0 + y ub1) = w0ᵀub0 + y(w0ᵀub1 + w1ᵀub0) + y²·w1ᵀub1 ; y² = −p·y − s
            t00 = (transpose(w0)*ub0)[1,1]; t01 = (transpose(w0)*ub1)[1,1]
            t10 = (transpose(w1)*ub0)[1,1]; t11 = (transpose(w1)*ub1)[1,1]
            C[1, j] = t00 - s*t11
            C[2, j] = t01 + t10 - p*t11
        end
        kc, KC = nullspace(C)
        for cc in 1:kc
            V = zero_matrix(Qe, 2n, 1)
            for j in 1:kv
                iszero(KC[j, cc]) && continue
                for i in 1:2n; V[i,1] += KC[j,cc]*KV[i,j]; end
            end
            iszero(V) && continue
            V = _qe_clear_denoms(ctx, V)
            w0 = V[1:n, 1:1]; w1 = V[(n+1):2n, 1:1]
            # s(y) = v(y)ᵀ u(y) = s0 + y s1 ; unit test via the field norm
            s0 = (transpose(w0)*u0)[1,1] - s*(transpose(w1)*u1)[1,1]
            s1 = (transpose(w0)*u1)[1,1] + (transpose(w1)*u0)[1,1] - p*(transpose(w1)*u1)[1,1]
            Nrm = s0^2 - p*s0*s1 + s*s1^2
            iszero(Nrm) && continue
            M0 = u0*transpose(w0) - s*(u1*transpose(w1))
            M1 = u0*transpose(w1) + u1*transpose(w0) - p*(u1*transpose(w1))
            Pi0 = (M0*(s0 - p*s1) + (s*s1)*M1) * inv(Nrm)
            Pi1 = (s0*M1 - s1*M0) * inv(Nrm)
            P = 2*Pi0 - p*Pi1
            (P*P == P) || continue
            rank(P) == 2 || continue
            push!(out, P)
        end
    end
    return out
end

# One Case A'' attempt.  Returns (moved, Anew, (T, desc)) like _one_balance!.
function _orbit_balance!(ctx::Ctx, A; verbose::Bool=false)
    n = size(A,1)
    factors = _alg_orbit_factors(ctx, A)
    isempty(factors) && return (false, A, nothing)
    obs0 = _obstruction_norm_alg(ctx, A)
    pts, _ = singular_points_x(ctx, A)
    polelist = collect(pts)
    # dump-point residues (full Q(ε))
    resid = Dict{Any,Any}()
    for p in polelist; resid[p] = residue_matrix(ctx, A, p); end
    resid[:inf] = residue_at_infinity(ctx, A)
    A2inf = nothing                                # lazy: 1/x² coefficient at ∞
    ncand = 0
    for cs in factors
        d = length(cs) - 1
        Rj = _orbit_residue_components(ctx, A, cs)
        Rj === nothing && continue
        all(iszero, Rj) && continue                # apparent orbit, no residue
        RjT = [transpose(R) for R in Rj]
        # rational eigenvalue candidates with a ≠ 0 (rational-eigvec branches satisfy
        # R₀u = λu, so λ must be an eigenvalue of R₀)
        for (lam, a, _) in _qe_linear_eigvals(ctx, Rj[1])
            iszero(a) && continue
            lower = a > 0                          # lower at roots (mirror) vs raise
            side_vecs = lower ? _orbit_joint_eigvecs(ctx, RjT, lam) :
                                _orbit_joint_eigvecs(ctx, Rj,  lam)
            isempty(side_vecs) && continue
            for D in _orbit_dump_configs(polelist, d)
                ncand > _ORBIT_MAXCAND && return (false, A, nothing)
                dpts = unique(D)
                infmult = count(p -> p === :inf, D)
                # simultaneous eigvec spaces at the dump points (left for raise, right
                # for lower — Lee's invariance condition at the shifted-down/up points)
                combos = Any[Any[]]                 # list of (point, μ) assignments
                okcfg = true
                for p in dpts
                    evs = _qe_linear_eigvals(ctx, resid[p])
                    isempty(evs) && (okcfg = false; break)
                    newcombos = Any[]
                    for cmb in combos, (mu, _, _) in evs
                        push!(newcombos, vcat(cmb, Any[(p, mu)]))
                    end
                    combos = newcombos
                end
                okcfg || continue
                for cmb in combos
                    ncand += 1
                    # stacked nullspace for the simultaneous dump eigvecs
                    M = nothing
                    for (p, mu) in cmb
                        R = resid[p]
                        B = lower ? (R - mu*identity_matrix(ctx.Qe, n)) :
                                    (transpose(R) - mu*identity_matrix(ctx.Qe, n))
                        M = M === nothing ? B : vcat(M, B)
                    end
                    k, K = nullspace(M)
                    k == 0 && continue
                    dump_vecs = Any[_qe_clear_denoms(ctx, K[:, j:j]) for j in 1:k]
                    # second-order condition at ∞ (∞ multiplicity 2 in D)
                    if infmult >= 2
                        if A2inf === nothing
                            A2inf = _laurent2_at_infinity(ctx, A)
                        end
                        W = transpose(dump_vecs[1])
                        for j in 2:length(dump_vecs); W = vcat(W, transpose(dump_vecs[j])); end
                        Mrefine = lower ? transpose(A2inf) : A2inf
                        refined = _refine_left_eig_in_rowspace(ctx, W, Mrefine)
                        dump_vecs = Any[transpose(v) for v in refined]
                        isempty(dump_vecs) && continue
                    end
                    for uvec in side_vecs, wvec in dump_vecs
                        # raise: u = orbit eigvec (column), v = dump eigvec; lower: swap
                        u  = lower ? wvec : uvec
                        vT = lower ? transpose(uvec) : transpose(wvec)
                        _projector_eps_regular(ctx, u, vT) || continue
                        P = (u * vT) * inv((vT*u)[1,1])
                        (P*P == P) || continue
                        T = _orbit_gauge(ctx, P, cs, D; lower_at_roots=lower)
                        Anew = try
                            apply_gauge(ctx, A, T)
                        catch
                            continue
                        end
                        is_fuchsian(ctx, Anew) || continue
                        _qex_mat_eps_regular(ctx, Anew) || continue
                        obs1 = _obstruction_norm_alg(ctx, Anew)
                        if obs1 < obs0
                            qstr = join([string(cs[k+1])*"*x^"*string(k) for k in 0:(d)], "+")
                            desc = "conjugate-orbit Lee balance: " *
                                   (lower ? "lower" : "raise") * " λ=$(lam) at the $(d) roots " *
                                   "of $(qstr), dump at $(join(string.(D), ","))"
                            verbose && println(stderr, "  [normalize] ", desc, "  obs(alg) ",
                                               obs0, " -> ", obs1)
                            return (true, Anew, (T, desc))
                        end
                    end
                end
            end
        end
        # --- Tr-pair (rank-2) variant for d = 2 (task-spec P = P_r + P_r̄): needed when the
        # orbit eigenvectors are genuinely K-valued (no joint rational eigvec, e.g. sec231's
        # x²+2x−14 with λ=−1−2ε).  λ candidates come from the hat model (full orbit spectrum).
        if d == 2
            Rhat = _orbit_hat_matrix(ctx, Rj, cs)
            for (lam, a, _) in _qe_linear_eigvals(ctx, Rhat)
                iszero(a) && continue
                lower = a > 0
                Ps = _orbit_trpair_projectors(ctx, Rj, cs, lam)
                isempty(Ps) && continue
                for D in _orbit_dump_configs(polelist, d), P in Ps
                    ncand += 1
                    ncand > _ORBIT_MAXCAND && return (false, A, nothing)
                    T = _orbit_gauge(ctx, P, cs, D; lower_at_roots=lower)
                    Anew = try
                        apply_gauge(ctx, A, T)
                    catch
                        continue
                    end
                    is_fuchsian(ctx, Anew) || continue
                    _qex_mat_eps_regular(ctx, Anew) || continue
                    obs1 = _obstruction_norm_alg(ctx, Anew)
                    if obs1 < obs0
                        qstr = join([string(cs[k+1])*"*x^"*string(k) for k in 0:d], "+")
                        desc = "conjugate-pair Tr-balance (rank-2): " *
                               (lower ? "lower" : "raise") * " λ=$(lam) at both roots " *
                               "of $(qstr), dump at $(join(string.(D), ","))"
                        verbose && println(stderr, "  [normalize] ", desc, "  obs(alg) ",
                                           obs0, " -> ", obs1)
                        return (true, Anew, (T, desc))
                    end
                end
            end
        end
        # --- Same-orbit ± shear (d = 2): the single rational orbit move T_q = M(x)/q̂. ---
        # When the orbit itself carries BOTH signs of the ε⁰ eigenvalue (a ± pair at the
        # conjugate roots), every move above transfers a-units between the orbit and OUTSIDE
        # dump points and cannot fire; per-root conjugate balances compose to an IRRATIONAL
        # gauge.  The rational cure is one matrix move
        #     T = I + N(x)/q̂,   N = (2x+p)·C₀ − (p·x+2s)·C₁,   C(y) = γ(y)·u(y)v(y)ᵀ,
        # i.e. T = M(x)/q̂ with M = q̂·I + N of x-degree ≤ 2 and det T = 1 (det M = q̂ⁿ):
        #   * u(y) is a K-eigvector of the orbit residue at the RAISED eigval λ_lo (a < 0),
        #     v(y)ᵀ a left K-eigvector at the LOWERED eigval λ_hi (a > 0), from the hat model;
        #     v(y)ᵀu(y) = 0 is automatic (biorthogonality at distinct eigvals over the field K);
        #   * cross-orthogonality v(y)ᵀu(ȳ) = 0 (imposed, as in the Tr-pair construction) kills
        #     the conjugate cross-products, so T is rational with T⁻¹ = 2I − T — equivalently
        #     adj(M) annihilates the balanced line to second order at both roots;
        #   * γ(y) = (λ_hi − λ_lo − 1) / (v(y)ᵀ B(y) u(y)) cancels the induced double pole at
        #     the roots (B = constant orbit Laurent term).  A nonzero-odd obstruction generally
        #     fails the cross gate — the algebraic-letter (irremovable) class stays a stop.
        # Result at both roots: λ_lo-line → λ_hi − 1, λ_hi-line → λ_lo + 1 — the same-point ±
        # transfer of Case A′, entirely inside the orbit, with no effect at ∞ or elsewhere.
        # Hard-gated like every other move: Fuchsian + ε-regular + strict drop of the
        # algebraic-aware obstruction norm.
        if d == 2
            p2 = ctx.Qe(cs[2]); s2 = ctx.Qe(cs[1])
            BK = _orbit_const_term(ctx, A, cs)
            Rhat  = _orbit_hat_matrix(ctx, Rj, cs)
            RhatT = _orbit_hat_matrix(ctx, RjT, cs)
            evs = BK === nothing ? Any[] : _qe_linear_eigvals(ctx, Rhat)
            for (lam_lo, a_lo, _) in evs, (lam_hi, a_hi, _) in evs
                (a_lo < 0 && a_hi > 0) || continue
                gnum = lam_hi - lam_lo - ctx.Qe(1)
                iszero(gnum) && continue               # γ = 0 would be the identity move
                ku, KU = nullspace(Rhat  - lam_lo*identity_matrix(ctx.Qe, 2n))
                kv, KV = nullspace(RhatT - lam_hi*identity_matrix(ctx.Qe, 2n))
                (ku == 0 || kv == 0) && continue
                for cu in 1:ku
                    U = _qe_clear_denoms(ctx, KU[:, cu:cu])
                    u0 = U[1:n, 1:1]; u1 = U[(n+1):2n, 1:1]
                    # conjugate ū(y) = (u₀ − p·u₁) − y·u₁; select v with v(y)ᵀ·ū(y) = 0 in K
                    ub0 = u0 - p2*u1; ub1 = -u1
                    Ccross = zero_matrix(ctx.Qe, 2, kv)
                    for jv in 1:kv
                        w0 = KV[1:n, jv:jv]; w1 = KV[(n+1):2n, jv:jv]
                        t00 = (transpose(w0)*ub0)[1,1]; t01 = (transpose(w0)*ub1)[1,1]
                        t10 = (transpose(w1)*ub0)[1,1]; t11 = (transpose(w1)*ub1)[1,1]
                        Ccross[1, jv] = t00 - s2*t11
                        Ccross[2, jv] = t01 + t10 - p2*t11
                    end
                    kc, KC = nullspace(Ccross)
                    for cc in 1:kc
                        ncand += 1
                        ncand > _ORBIT_MAXCAND && return (false, A, nothing)
                        V = zero_matrix(ctx.Qe, 2n, 1)
                        for jv in 1:kv
                            iszero(KC[jv, cc]) && continue
                            for iv in 1:2n; V[iv,1] += KC[jv,cc]*KV[iv,jv]; end
                        end
                        iszero(V) && continue
                        V = _qe_clear_denoms(ctx, V)
                        v0 = V[1:n, 1:1]; v1 = V[(n+1):2n, 1:1]
                        # b = v(y)ᵀ B(y) u(y) ∈ K  (K-trilinear contraction, reduced mod q̂)
                        B0, B1 = BK
                        Bu0 = B0*u0 - s2*(B1*u1)
                        Bu1 = B0*u1 + B1*u0 - p2*(B1*u1)
                        b0 = (transpose(v0)*Bu0)[1,1] - s2*(transpose(v1)*Bu1)[1,1]
                        b1 = (transpose(v0)*Bu1)[1,1] + (transpose(v1)*Bu0)[1,1] -
                             p2*(transpose(v1)*Bu1)[1,1]
                        binv = _k2_inv(p2, s2, (b0, b1))
                        binv === nothing && continue    # double-pole condition unsolvable
                        g0 = gnum*binv[1]; g1 = gnum*binv[2]      # γ(y) = γ₀ + y·γ₁
                        # C(y) = γ(y)·u(y)v(y)ᵀ = C₀ + y·C₁ (invariant under K-rescaling of u, v)
                        W0 = u0*transpose(v0) - s2*(u1*transpose(v1))
                        W1 = u0*transpose(v1) + u1*transpose(v0) - p2*(u1*transpose(v1))
                        C0 = g0*W0 - s2*(g1*W1)
                        C1 = g0*W1 + g1*W0 - p2*(g1*W1)
                        # T = I + [(2x+p)·C₀ − (p·x+2s)·C₁] / q̂
                        x = ctx.x
                        qx = sum(ctx.Qex(cs[k+1]) * x^k for k in 0:d)
                        Nmat = (2*x + ctx.Qex(p2)) * Fuchsia._qe_to_qex(ctx, C0) -
                               (ctx.Qex(p2)*x + 2*ctx.Qex(s2)) * Fuchsia._qe_to_qex(ctx, C1)
                        T = identity_matrix(ctx.Qex, n) + (ctx.Qex(1)//qx) * Nmat
                        Anew = try
                            apply_gauge(ctx, A, T)
                        catch
                            continue
                        end
                        is_fuchsian(ctx, Anew) || continue
                        _qex_mat_eps_regular(ctx, Anew) || continue
                        obs1 = _obstruction_norm_alg(ctx, Anew)
                        if obs1 < obs0
                            qstr = join([string(cs[k+1])*"*x^"*string(k) for k in 0:d], "+")
                            desc = "same-orbit ± shear T_q=M(x)/q: raise λ=$(lam_lo), " *
                                   "lower λ=$(lam_hi) at both roots of $(qstr)"
                            verbose && println(stderr, "  [normalize] ", desc, "  obs(alg) ",
                                               obs0, " -> ", obs1)
                            return (true, Anew, (T, desc))
                        end
                    end
                end
            end
        end
    end
    return (false, A, nothing)
end

function _constant_epsdecouple(ctx::Ctx, A)
    # extract eps^0 connection A0(x)
    _, comps = eps_decompose(ctx, A)
    A0 = comps[1]
    if iszero(A0)
        return (false, nothing, A)   # nothing to do here
    end
    # Try to remove A0 by a rational gauge W with dW = A0 W, W normalized at infinity = I.
    Wok, W = _rational_coboundary(ctx, A0)
    Wok || return (false, nothing, A)
    Winv = inv(W)
    # apply T = Winv as a gauge to the FULL A (includes the dT term)
    Anew = apply_gauge(ctx, A, Winv)
    return (true, Winv, Anew)
end

# Solve dW/dx = A0 W for a rational matrix W (entries rational in x), W(∞) = I, when A0 has
# nilpotent residues so that W is rational (no genuine logs).  Method: Frobenius/recursion in
# the partial-fraction basis.  We expand W as a finite sum over products of the simple-pole
# generators g_k = 1/(x - x_k) up to total degree L (nilpotency bound), solving for constant
# matrix coefficients by matching dW = A0 W coefficient-wise.  Returns (ok, W).
function _rational_coboundary(ctx::Ctx, A0; maxL::Int=8)
    n = size(A0,1)
    pts, _ = singular_points_x(ctx, A0)
    poles = collect(pts)
    isempty(poles) && return (true, identity_matrix(ctx.Qex, n))
    # residue matrices N_k (constant) of A0
    Nk = Dict(p => eps0_part(ctx, residue_matrix(ctx, A0, p)) for p in poles)
    # Strict triangularizability test: all N_k must be nilpotent (needed for rational W).
    for (p, N) in Nk
        if !_is_nilpotent(N)
            return (false, nothing)   # genuine eps^0 eigenvalue -> not a coboundary; balances handle it
        end
    end
    # Ansatz: W(x) = I + Σ over nonempty words w=(k_1,...,k_r), r<=L, of  M_w * f_w(x),
    # where f_w(x) are the iterated-pole functions defined by  f_∅ = 1,  d/dx f_{k·w'} = g_k f_{w'}
    # with g_k = 1/(x-x_k).  Then  dW = Σ M_w (Σ_k? )...  Equivalently we solve the recursion
    #   d/dx W = A0 W  =>  matching: the coefficient structure is the path-ordered (Chen) series,
    #   which TERMINATES because the N_k are nilpotent and (for a coboundary) the higher Chen
    #   integrals of nilpotent letters vanish in the rational closure.
    # We realize it by direct truncated Picard iteration in the rational function field:
    #     W_0 = I;  W_{m+1} = I + ∫ A0 W_m dx   (∫ taken in the rational closure, i.e. only the
    #   exact-rational antiderivative; if a genuine log appears, FAIL -> not rational coboundary).
    W = identity_matrix(ctx.Qex, n)
    for _ in 1:maxL
        rhs = A0 * W                      # = dW_target
        intg, ok = _rational_integral_matrix(ctx, rhs, poles)
        ok || return (false, nothing)
        Wnew = identity_matrix(ctx.Qex, n) + intg
        if Wnew == W
            # converged: verify dW = A0 W exactly
            if dmat_dx(ctx, W) == A0 * W
                return (true, W)
            else
                return (false, nothing)
            end
        end
        W = Wnew
    end
    # final verify after maxL iterations
    if dmat_dx(ctx, W) == A0 * W
        return (true, W)
    end
    return (false, nothing)
end

# nilpotency test for a constant QQ matrix
function _is_nilpotent(N)
    cp = charpoly(N)
    n = size(N,1)
    # nilpotent <=> charpoly = x^n
    return cp == gen(parent(cp))^n
end

# Rational antiderivative of a Qex matrix whose entries are sums of c/(x-x_k) (+ polynomial),
# vanishing at infinity (no constant of integration; we want the piece that ->0 at infinity).
# Returns (integral_matrix, ok).  ok=false if a 1/(x-x_k) term is present (would give a log).
function _rational_integral_matrix(ctx::Ctx, M, poles)
    m, n = size(M)
    out = zero_matrix(ctx.Qex, m, n)
    for i in 1:m, j in 1:n
        f = M[i,j]
        iszero(f) && continue
        g, ok = _rational_integral_scalar(ctx, f, poles)
        ok || return (out, false)
        out[i,j] = g
    end
    return (out, true)
end

# antiderivative of a scalar Qex rational function, requiring NO simple-pole (log) residues.
# Uses partial fractions: f = poly(x) + Σ_k [ a_{k,1}/(x-x_k) + a_{k,2}/(x-x_k)^2 + ... ].
# The 1/(x-x_k) coefficients must vanish (else log); higher powers integrate to rationals;
# the polynomial part integrates to a polynomial.  We also kill the polynomial part's constant
# of integration by demanding the result -> 0 at infinity (drop the constant).
function _rational_integral_scalar(ctx::Ctx, f, poles)
    x = ctx.x
    num = numerator(f); den = denominator(f)
    # residue at each pole must be zero (no log)
    for p in poles
        r = _scalar_residue(ctx, f, p)
        iszero(r) || return (ctx.Qex(0), false)
    end
    # also residue at infinity must be zero for a rational antiderivative vanishing at infinity
    # (the polynomial part is allowed; it gives a higher polynomial).  Build the antiderivative
    # by integrating the partial-fraction decomposition with zero log residues.
    # Easiest exact route: since all simple-pole residues vanish, f has only poles of order>=2
    # plus a polynomial.  Its antiderivative is rational.  Compute via the Hermite reduction:
    #   antiderivative of 1/(x-a)^m (m>=2) = -1/((m-1)(x-a)^{m-1}); of x^j = x^{j+1}/(j+1).
    g, ok = _hermite_antideriv(ctx, f, poles)
    return (g, ok)
end

function _scalar_residue(ctx::Ctx, f, p)
    # reuse RatFunc residue machinery
    return RatFunc._residue_scalar(ctx, f, ctx.Qe(p))
end

# Hermite-style exact antiderivative for a rational function whose only obstruction (log) terms
# (simple-pole residues) already vanish.  We compute it by a partial-fraction expansion over the
# given poles + the polynomial quotient.
function _hermite_antideriv(ctx::Ctx, f, poles)
    x = ctx.x
    num = numerator(f); den = denominator(f)
    Px = parent(num)
    xx = gen(Px)
    # polynomial part q and remainder r/den
    q, r = divrem(num, den)
    # integrate polynomial part q -> Q (drop constant)
    Q = ctx.Qex(0)
    for i in 0:degree(q)
        ci = coeff(q, i)
        Q += ctx.Qex(ci) * x^(i+1) // ctx.Qex(i+1)
    end
    # partial fractions of r/den over the poles: for each pole p of multiplicity μ, get
    # coefficients a_{p,1..μ}; a_{p,1} (the log term) is zero by assumption.
    acc = Q
    rem_f = ctx.Qex(r) // ctx.Qex(den)
    for p in poles
        μ = _pole_mult(ctx, rem_f, p)
        μ == 0 && continue
        # extract a_{p,m} for m=1..μ via Laurent coefficients: a_{p,m} = res of (x-p)^{m-1} f?
        # use: a_{p,m} = (1/(μ-m)!) d^{μ-m}/dx^{μ-m} [ (x-p)^μ f ] |_{x=p}
        for m in 1:μ
            a = _laurent_coeff(ctx, rem_f, p, m, μ)
            iszero(a) && continue
            if m == 1
                # log term — must be zero
                iszero(a) || return (ctx.Qex(0), false)
            else
                # ∫ a/(x-p)^m = -a/((m-1)(x-p)^{m-1})
                acc += -ctx.Qex(a) // (ctx.Qex(m-1) * (x - ctx.Qex(p))^(m-1))
            end
        end
    end
    return (acc, true)
end

function _pole_mult(ctx::Ctx, f, p)
    den = denominator(f)
    Px = parent(den)
    xx = gen(Px)
    lin = xx - ctx.Qe(p)
    k = 0; dd = den
    while true
        qd, rd = divrem(dd, lin)
        if iszero(rd); dd = qd; k += 1 else break end
    end
    return k
end

# Laurent coefficient a_{p,m} (coeff of (x-p)^{-m}) of a rational f with pole order μ at p.
function _laurent_coeff(ctx::Ctx, f, p, m, μ)
    # a_{p,m} = 1/(μ-m)! * d^{μ-m}/dx^{μ-m} [ (x-p)^μ f ] |_{x=p}
    x = ctx.x
    h = ((x - ctx.Qex(p))^μ) * f
    for _ in 1:(μ-m)
        h = derivative(h)
    end
    val = RatFunc._eval_qex_at(ctx, h, ctx.Qe(p))
    return val // ctx.Qe(factorial(big(μ-m)))
end

# --------------------------------------------------------------------------------------------
# A single eps^0-reducing balance.  Strategy (Lee, semisimple case):
#   * Look at the residue at each singular point.  Its eps^0 part R0 = res|_{eps=0}.
#   * If some point x1 has R0 with a nonzero eigenvalue, balance it against another point x2
#     (preferring infinity) whose R0 has a compensating eigenvalue, using the spectral
#     projector onto that eigenvalue.  After the balance the offending eps^0 eigenvalues are
#     shifted toward zero.
# We implement the common special case arising in UT systems: the eps^0 part is a
# constant nilpotent/decoupling block removable by a balance with a projector built from the
# generalized eigenspace.  We construct projectors over Q (eps=0 specialization) and lift.
# --------------------------------------------------------------------------------------------

function _one_balance!(ctx::Ctx, A, A0; verbose::Bool=false)
    n = size(A,1)
    pts, _ = singular_points_x(ctx, A)
    allpts = vcat(collect(pts), [:inf])
    # FULL Q(eps) residues (the rank-1 projector must be built from Q(eps) eigenvectors so that
    # the balance preserves Fuchsianity exactly — using eps=0 eigenvectors would only kill the
    # double-pole term to O(eps^0), leaving an O(eps) order-2 pole).
    resid = Dict{Any,Any}()
    for a in pts; resid[a] = residue_matrix(ctx, A, a); end
    resid[:inf] = residue_at_infinity(ctx, A)
    res0 = Dict(p => eps0_part(ctx, resid[p]) for p in allpts)

    # ----------------------------------------------------------------------------------------
    # Lee normalization (1411.0911 §4), CORRECTED rank-1 balance.
    #
    # For B(P, x1, x2) = (I-P) + (x-x2)/(x-x1) P with rank-1  P = u vᵀ / (vᵀ u), the balance
    # preserves Fuchsianity iff
    #     P R(x1) (I-P) = 0   ⇔  vᵀ is a LEFT  eigenvector of the FULL residue R(x1),
    #     (I-P) R(x2) P = 0   ⇔  u  is a RIGHT eigenvector of the FULL residue R(x2),
    # and then the eigenvalue carried by vᵀ at x1 is LOWERED by 1 and the eigenvalue carried by
    # u at x2 is RAISED by 1.  We therefore:
    #   * collect every (pole, λ(eps), a = eps^0 part) eigenvalue over ALL singular points,
    #   * pick x1=posp with a>0 (left eigvec vᵀ) and x2=negp with a<0 (right eigvec u), posp≠negp,
    #   * build P over Q(eps) and balance — guaranteed Σ|a| drop of 2 per round.
    #
    # The previous implementation (a) greedily set posp=negp = the FIRST pole that happened to
    # carry both signs, then skipped Case A; and (b) used the full SPECTRAL projector at one
    # point, which generically violates the second invariance condition.  Both are fixed here.
    # ----------------------------------------------------------------------------------------

    pos = Any[]   # (pole, λ∈Qe, a∈QQ, vᵀ row-vec over Qe)   — to be LOWERED at x1
    neg = Any[]   # (pole, λ∈Qe, a∈QQ, u  col-vec over Qe)   — to be RAISED  at x2
    zro = Any[]   # (pole, λ∈Qe,        u, vᵀ)               — a=0 partners for same-pole routing
    for p in allpts
        R = resid[p]
        for (lam, a, e) in _qe_linear_eigvals(ctx, R)
            if a > 0
                for vT in _qe_left_eigvecs(ctx, R, lam)
                    push!(pos, (p, lam, a, vT))
                end
            elseif a < 0
                for u in _qe_right_eigvecs(ctx, R, lam)
                    push!(neg, (p, lam, a, u))
                end
            else
                # keep one a=0 direction per pole for same-pole routing fallback
                us = _qe_right_eigvecs(ctx, R, lam)
                vs = _qe_left_eigvecs(ctx, R, lam)
                if !isempty(us) && !isempty(vs)
                    push!(zro, (p, lam, us[1], vs[1]))
                end
            end
        end
    end
    # greedy on the largest off-target integer part first
    sort!(pos; by = t -> -t[3])
    sort!(neg; by = t ->  t[3])

    obs0 = _obstruction_norm(ctx, A)
    _pt(p) = p === :inf ? :inf : ctx.Qe(p)

    # ---- Case A (primary): cross-pole Lee rank-1 balance. ----
    for (pp, lpos, ap, vT) in pos, (np, lneg, an, u) in neg
        pp === np && continue
        _projector_eps_regular(ctx, u, vT) || continue
        P = (u * vT) * inv((vT*u)[1,1])
        T = gauge_balance(ctx, P, _pt(pp), _pt(np))
        Anew = apply_gauge(ctx, A, T)
        is_fuchsian(ctx, Anew) || continue
        if _obstruction_norm(ctx, Anew) < obs0
            desc = "Lee rank-1 balance: lower λ=$(lpos)@$(pp) raise λ=$(lneg)@$(np)"
            verbose && println(stderr, "  [normalize] ", desc, "  obs ", obs0, " -> ",
                               _obstruction_norm(ctx, Anew))
            return (true, Anew, (T, desc))
        end
    end

    # ---- Case A' (same-pole ±, routed): posp==negp via an intermediate partner q. ----
    # Composite of two rank-1 balances B(P2, q, p) ∘ B(P1, p, q) lowers a⁺@p AND raises a⁻@p
    # while leaving q's eps^0 spectrum net-unchanged.  Needed when every nonzero a is at one pole.
    # The two projectors are GENERICALLY eps-singular (eigvecs at q for a=0 eigvals collapse at
    # eps→0), but the COMPOSITE gauge T2*T1 produces an eps-regular connection — the 1/eps poles
    # of P1 and P2 cancel.  We therefore drop the eps-regularity gate per-step and instead check
    # eps-regularity + Fuchsianity + obstruction-drop on the composite result A2 only.
    for (pp, lpos, ap, vT) in pos, (np, lneg, an, _) in neg
        pp === np || continue
        for (qp, lq, uq, _) in zro
            qp === pp && continue
            s1 = (vT*uq)[1,1]
            iszero(s1) && continue
            P1 = (uq * vT) * inv(s1)
            T1 = gauge_balance(ctx, P1, _pt(pp), _pt(qp))
            A1 = apply_gauge(ctx, A, T1)
            is_fuchsian(ctx, A1) || continue
            # step 2: re-derive eigvecs from A1's residues (step-1 is NOT a similarity at p,q).
            Rq1 = (qp === :inf) ? residue_at_infinity(ctx, A1) : residue_matrix(ctx, A1, qp)
            Rp1 = (np === :inf) ? residue_at_infinity(ctx, A1) : residue_matrix(ctx, A1, np)
            for (lq2, aq2, _) in _qe_linear_eigvals(ctx, Rq1)
                aq2 > 0 || continue
                for vT2 in _qe_left_eigvecs(ctx, Rq1, lq2),
                    (ln2, an2, _) in _qe_linear_eigvals(ctx, Rp1)
                    an2 < 0 || continue
                    for u2 in _qe_right_eigvecs(ctx, Rp1, ln2)
                        s2 = (vT2*u2)[1,1]
                        iszero(s2) && continue
                        P2 = (u2 * vT2) * inv(s2)
                        T2 = gauge_balance(ctx, P2, _pt(qp), _pt(np))
                        A2 = apply_gauge(ctx, A1, T2)
                        is_fuchsian(ctx, A2) || continue
                        _qex_mat_eps_regular(ctx, A2) || continue   # composite must be eps-regular
                        if _obstruction_norm(ctx, A2) < obs0
                            Tcomp = T2 * T1
                            desc = "same-pole ± routed via $(qp): lower λ=$(lpos)@$(pp), " *
                                   "raise λ=$(lneg)@$(np)"
                            verbose && println(stderr, "  [normalize] ", desc, "  obs ", obs0,
                                               " -> ", _obstruction_norm(ctx, A2))
                            return (true, A2, (Tcomp, desc))
                        end
                    end
                end
            end
        end
    end

    # ---- Legacy fallback (kept intact): single-point spectral-projector balance. ----
    posp = nothing; pos_lam = nothing; pos_P = nothing
    negp = nothing; neg_lam = nothing
    for p in allpts
        M = res0[p]; iszero(M) && continue
        for r in rational_roots(charpoly_q(M))
            if r > 0 && posp === nothing
                P = spectral_projector(M, r)
                P !== nothing && (posp = p; pos_lam = r; pos_P = P)
            elseif r < 0 && negp === nothing
                negp = p; neg_lam = r
            end
        end
    end
    if posp !== nothing
        for q in allpts
            q === posp && continue
            T = gauge_balance(ctx, pos_P, _pt(posp), _pt(q))
            Anew = apply_gauge(ctx, A, T)
            if is_fuchsian(ctx, Anew) && _obstruction_norm(ctx, Anew) < obs0
                return (true, Anew, (T, "legacy spectral balance lower@$(posp) raise@$(q)"))
            end
        end
    end
    if negp !== nothing
        Pneg = spectral_projector(res0[negp], neg_lam)
        if Pneg !== nothing
            for q in allpts
                q === negp && continue
                T = gauge_balance(ctx, Pneg, _pt(q), _pt(negp))
                Anew = apply_gauge(ctx, A, T)
                if is_fuchsian(ctx, Anew) && _obstruction_norm(ctx, Anew) < obs0
                    return (true, Anew, (T, "legacy spectral balance lower@$(q) raise@$(negp)"))
                end
            end
        end
    end

    # ---- Case A'' (conjugate-pair / Galois-orbit): balances at irreducible algebraic
    # letters (quadratic/cubic ε-free denominator factors).  Reached only when every
    # rational-point move above has failed, so pre-existing behaviour is unchanged. ----
    movedq, Aq, Tq = _orbit_balance!(ctx, A; verbose=verbose)
    movedq && return (true, Aq, Tq)

    return (false, A, nothing)
end

# --- Q(eps) eigen helpers (linear-factor route; Nemo factors over Q(eps)) -------------------

"Eigenvalues of a Qe matrix R as Q(eps) elements, via linear factors of charpoly.  Returns a
 list of (λ::Qe, a::QQ = eps^0 part of λ, mult::Int).  Non-linear irreducible factors (eigvals
 not in Q(eps)) are skipped — for Feynman/dlog residues every eigval is affine in eps."
function _qe_linear_eigvals(ctx::Ctx, R)
    Qe = ctx.Qe
    cp = charpoly(R)
    out = Tuple{Any,QQFieldElem,Int}[]
    for (q, e) in factor(cp)
        degree(q) == 1 || continue
        lam = -coeff(q,0)//coeff(q,1)
        a = Fuchsia._eval_qe_at0(lam)
        push!(out, (lam, a, e))
    end
    return out
end

"Right eigenvectors (columns over Qe) of R for eigval λ, with eps-denominators cleared."
function _qe_right_eigvecs(ctx::Ctx, R, lam)
    n = size(R,1)
    k, K = nullspace(R - lam*identity_matrix(ctx.Qe, n))
    return [_qe_clear_denoms(ctx, K[:, j:j]) for j in 1:k]
end

"Left eigenvectors (rows over Qe) of R for eigval λ, with eps-denominators cleared."
function _qe_left_eigvecs(ctx::Ctx, R, lam)
    n = size(R,1)
    k, K = nullspace(transpose(R) - lam*identity_matrix(ctx.Qe, n))
    return [_qe_clear_denoms(ctx, transpose(K[:, j:j])) for j in 1:k]
end

# Clear common eps-denominators (and any overall eps^k factor in the numerators) so the vector
# has polynomial-in-eps entries with at least one entry nonzero at eps=0.  Makes the projector
# P = u vᵀ/(vᵀu) well-defined and lets us test eps-regularity of P via vᵀu|_{eps=0}.
function _qe_clear_denoms(ctx::Ctx, V)
    m,n = size(V)
    Pe = parent(numerator(ctx.eps))     # Q[eps] poly ring
    L = one(Pe)
    for i in 1:m, j in 1:n
        iszero(V[i,j]) && continue
        L = lcm(L, Pe(denominator(V[i,j])))
    end
    W = ctx.Qe(L) * V
    # strip any overall common eps-polynomial factor in numerators (so some entry has eps^0 ≠ 0)
    g = zero(Pe)
    for i in 1:m, j in 1:n
        iszero(W[i,j]) && continue
        g = iszero(g) ? Pe(numerator(W[i,j])) : gcd(g, Pe(numerator(W[i,j])))
    end
    if !iszero(g) && degree(g) >= 1
        W = W * inv(ctx.Qe(g))
    end
    return W
end

# Is a Qex matrix regular at eps=0 (no eps in any x-coefficient denominator)?
function _qex_mat_eps_regular(ctx::Ctx, A)
    m,n = size(A)
    for i in 1:m, j in 1:n
        a = A[i,j]; iszero(a) && continue
        for p in (numerator(a), denominator(a))
            for k in 0:degree(p)
                c = coeff(p,k); iszero(c) && continue
                iszero(coeff(denominator(c),0)) && return false
            end
        end
    end
    return true
end

# eps^0 of a Qe element, or `nothing` if eps=0 is a pole.
function _safe_eps0(v)
    den = denominator(v)
    iszero(coeff(den,0)) && return nothing
    return coeff(numerator(v),0) // coeff(den,0)
end

# is the rank-1 projector P = u vᵀ/(vᵀu) regular and nonzero at eps=0?
function _projector_eps_regular(ctx::Ctx, u, vT)
    s = (vT * u)[1,1]
    iszero(s) && return false
    s0 = _safe_eps0(s)
    s0 === nothing && return false       # vᵀu has an eps pole (shouldn't after clearing)
    if iszero(s0)
        # vᵀu vanishes at eps=0 ⇒ P has a 1/eps pole UNLESS every (u_i v_j)|_{eps=0}=0 too.
        # Clearing denoms guarantees u,v each have a nonzero eps^0 entry, so generically P is
        # eps-singular here.  Skip; another (u,v) pair will work.
        return false
    end
    return true
end

# obstruction measure = sum over all singular points (incl. infinity) of the sum of |off-target
# eps^0 eigenvalue| of the residue.  Lee's normalization strictly decreases this.  We also add
# the count of nonzero eps^0 entries as a tiebreaker so that a final constant rotation (which
# zeroes the eps^0 part without changing eigenvalues) is recognized as progress.
function _obstruction_norm(ctx::Ctx, A)
    pts, _ = singular_points_x(ctx, A)
    tot = QQ(0)
    for p in vcat(collect(pts), [:inf])
        R = (p === :inf) ? residue_at_infinity(ctx, A) : residue_matrix(ctx, A, p)
        R0 = eps0_part(ctx, R)
        cp = charpoly_q(R0)
        for r in rational_roots_with_mult(cp)
            tot += abs(r)
        end
    end
    # entry count tiebreaker, scaled small
    _, comps = eps_decompose(ctx, A)
    A0 = comps[1]
    cnt = 0
    m, n = size(A0)
    for i in 1:m, j in 1:n
        iszero(A0[i,j]) || (cnt += 1)
    end
    return tot + QQ(cnt, 1000000)
end

# rational roots with multiplicity
function rational_roots_with_mult(p)
    rts = QQFieldElem[]
    fp = factor(p)
    for (q, e) in fp
        if degree(q) == 1
            c1 = coeff(q,1); c0 = coeff(q,0)
            for _ in 1:e
                push!(rts, -c0//c1)
            end
        end
    end
    return rts
end

# rank of a rational matrix
function rank_q(P)
    return rank(P)
end

# Build a projector from a constant (Qe but eps^0 already) matrix M (Q-rational here).
# Returns (P, eigenvalue) where P projects onto an eigenspace with NONZERO eigenvalue, or onto
# the image (nilpotent case).  Projector is a Q-rational matrix.
function _projector_from(M)
    n = size(M,1)
    # specialize: M is already eps^0 part, entries are Q-rationals (QQFieldElem)
    # eigen-decomposition over Q is generally not available; we use the spectral approach via
    # the minimal-poly free part.  Simplest robust move for the semisimple case:
    #   pick a nonzero column space vector and its dual to form a rank-1 projector along the
    #   dominant eigen direction.  For a clean balance we want P with M P = lam P + ... .
    # We implement: if M is diagonalizable over Q with a rational nonzero eigenvalue, build the
    # spectral projector; else use image/kernel complementary projector (nilpotent shift).
    QQm = M
    # try rational eigenvalues = rational roots of charpoly
    cp = charpoly_q(QQm)
    roots = rational_roots(cp)
    for r in roots
        iszero(r) && continue
        P = spectral_projector(QQm, r)
        P === nothing && continue
        return (P, r)
    end
    # nilpotent / no nonzero rational eigenvalue: project onto image of M (rank-deficient shift)
    P = image_projector(QQm)
    if P !== nothing && !iszero(P)
        return (P, QQ(0))
    end
    return (nothing, nothing)
end

# charpoly over QQ
function charpoly_q(M)
    return charpoly(M)
end

# rational roots of a QQ-univariate polynomial
function rational_roots(p)
    rts = QQFieldElem[]
    fp = factor(p)
    for (q, _) in fp
        if degree(q) == 1
            c1 = coeff(q,1); c0 = coeff(q,0)
            push!(rts, -c0//c1)
        end
    end
    return rts
end

# spectral projector onto eigenspace of eigenvalue r of M (over QQ), if M is semisimple on it:
#   P = product over other eigenvalues s of (M - s I)/(r - s), evaluated; valid when M acts
#   semisimply (we verify P^2 = P).
function spectral_projector(M, r)
    n = size(M,1)
    cp = charpoly_q(M)
    rts = rational_roots(cp)
    others = [s for s in unique(rts) if s != r]
    P = identity_matrix(QQ, n)
    for s in others
        P = P * (M - s*identity_matrix(QQ,n))
        denom = (r - s)
        P = P * inv(denom)
    end
    # normalize so that P is idempotent if M semisimple: with the full product over distinct
    # roots this is the Lagrange interpolation projector already.  Verify.
    if P*P == P
        return P
    end
    return nothing
end

# projector onto image(M) along a complement (for nilpotent shifts): returns M*pinv-like.
function image_projector(M)
    n = size(M,1)
    # columns of M span the image; build projector via RREF basis. Use M itself scaled if M^2 ~ M.
    # Simplest: if M is already idempotent-up-to-scale, normalize.  Otherwise return nothing
    # (defer to the STOP path — non-semisimple obstruction).
    if iszero(M); return nothing; end
    # check if M is a scalar multiple of an idempotent
    # find a nonzero entry to guess scale via trace of a rank-1: skip — return nothing to STOP.
    return nothing
end

end # module

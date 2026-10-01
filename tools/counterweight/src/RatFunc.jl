# RatFunc.jl — rational-function & matrix utilities for the eps-factorization engine.
#
# Representation (the key design choice):
#   * eps lives in a base rational-function field  Qe = Q(eps).
#   * the distinguished differentiation variable x lives in  Qex = Qe(x).
#   So a connection matrix A(eps,x) is a Nemo matrix whose entries are univariate rational
#   functions in x with coefficients in Q(eps).  This makes pole/residue analysis in x exact
#   and trivial (factor the denominator), which is what Lee's algorithm needs.
#
# Any further kinematic variables (e.g. a second ratio y) are handled by Counterweight.jl:
#   the eps-form is constructed in x with y frozen at exact rationals, then the SAME constant
#   gauge T is checked to eps-factorize the y-connection too (integrability / flatness).
#
# No Wolfram, no Sage.  Pure Julia + Nemo (Flint).

module RatFunc

using Nemo

export Ctx, make_ctx, F2x, matF2x,
       singular_points_x, residue_matrix, residue_at_infinity,
       is_fuchsian, poly_pole_order, dmat_dx, eval_x_at,
       eps_of, x_of,
       pole_profile, laurent_coeff_matrix, subst_x_inv, subst_x_inv_mat

"""
    Ctx

Field tower for the engine:
  Qe  = Q(eps)             (base, eps generator `eps`)
  Qex = Qe(x)              (x generator `x`)
Differentiation is d/dx on Qex.
"""
struct Ctx
    Qe
    eps
    Qex
    x
end

function make_ctx(; epsname::String="eps", xname::String="x")
    Qe, eps = rational_function_field(QQ, epsname)
    Qex, x = rational_function_field(Qe, xname)
    return Ctx(Qe, eps, Qex, x)
end

eps_of(ctx::Ctx) = ctx.eps
x_of(ctx::Ctx)   = ctx.x

"Coerce a number / Qe element / Qex element into Qex."
F2x(ctx::Ctx, v) = ctx.Qex(v)

"Build an n×n Qex matrix from a function entry(i,j) returning Qex-coercible values."
function matF2x(ctx::Ctx, n::Int, entry)
    A = zero_matrix(ctx.Qex, n, n)
    for i in 1:n, j in 1:n
        A[i,j] = F2x(ctx, entry(i,j))
    end
    return A
end

"""
    dmat_dx(ctx, A)

Entrywise d/dx of a Qex matrix.  Nemo provides `derivative` on rational_function_field elems.
"""
function dmat_dx(ctx::Ctx, A)
    m, n = size(A)
    B = zero_matrix(ctx.Qex, m, n)
    for i in 1:m, j in 1:n
        B[i,j] = derivative(A[i,j])
    end
    return B
end

"""
    singular_points_x(ctx, A) -> Vector of Qe-elements (finite poles in x)

Collect the finite singular points: roots of all entry denominators that are *linear* in x
with coefficients in Q(eps).  For Fuchsian dlog systems every letter is linear in x after the
ratio change of variables, so the singularities are rational points a(eps) ∈ Q(eps) — usually
plain rationals (0, 1, ...).  Returns the distinct pole locations as Qe elements.

Non-linear irreducible denominator factors (algebraic letters / square roots) are reported via
`algebraic_factors` (caller decides: rationalize first, or treat as a genuine algebraic letter).
"""
function singular_points_x(ctx::Ctx, A)
    x = ctx.x
    pts = Dict{Any,Bool}()   # use string key for dedup
    algfac = Any[]
    m, n = size(A)
    for i in 1:m, j in 1:n
        a = A[i,j]
        iszero(a) && continue
        d = denominator(a)
        is_constant(d) && continue
        fd = factor(d)
        for (p, _) in fd
            dp = degree(p)
            if dp == 1
                # p = c1*x + c0  -> root = -c0/c1 ∈ Qe
                c1 = coeff(p, 1); c0 = coeff(p, 0)
                root = -c0 // c1
                pts[string(root)] = true
                pts[string(root)] = true
                _store_point!(pts, root)
            elseif dp >= 2
                push!(algfac, p)
            end
        end
    end
    # rebuild list of unique Qe roots
    return _collect_points(pts), algfac
end

# helpers to dedup Qe roots by string but keep the value
const _PTVAL = Dict{String,Any}()
function _store_point!(pts, root)
    _PTVAL[string(root)] = root
    pts[string(root)] = true
end
function _collect_points(pts)
    out = Any[]
    for k in keys(pts)
        push!(out, _PTVAL[k])
    end
    return out
end

"""
    residue_matrix(ctx, A, a) -> matrix over Qe

Matrix residue of A(eps,x) at the simple/multiple finite pole x = a (a ∈ Qe).
Residue = coefficient of 1/(x-a) in the Laurent expansion, returned as a Qe-matrix.

Implementation: residue_{x=a} f = lim_{x->a} (x-a) f for a simple pole.  For safety against
higher-order spurious poles we extract via the Laurent series: shift x = a + t, take the
coefficient of t^{-1}.  We do this by computing ((x-a)^p * f) evaluated at x=a after enough
derivatives — but the clean way with Nemo is the partial-fraction-free residue:
    res = numerator(g)(a) / (d/dx denominator(g))(a)   when den has a simple root at a,
which we generalize by clearing the (x-a) multiplicity.
Here we use the robust Laurent route via Taylor coefficients.
"""
function residue_matrix(ctx::Ctx, A, a)
    m, n = size(A)
    Qe = ctx.Qe
    R = zero_matrix(Qe, m, n)
    for i in 1:m, j in 1:n
        R[i,j] = _residue_scalar(ctx, A[i,j], a)
    end
    return R
end

# residue of a scalar Qex element at x=a (a in Qe)
function _residue_scalar(ctx::Ctx, f, a)
    Qe = ctx.Qe; x = ctx.x
    iszero(f) && return Qe(0)
    num = numerator(f); den = denominator(f)
    # multiplicity of (x-a) in den
    p = x - ctx.Qex(a)   # but need polynomial form; work in poly ring
    # Move to the underlying univariate polynomial ring over Qe.
    Px = parent(num)              # poly ring Qe[x]
    aa = Qe(a)
    # factor out (x - aa)^k from den
    k = 0
    dd = den
    xx = gen(Px)
    lin = xx - aa
    while true
        q, r = divrem(dd, lin)
        if iszero(r)
            dd = q; k += 1
        else
            break
        end
    end
    if k == 0
        return Qe(0)   # no pole at a
    end
    # f = num / ((x-a)^k * dd).  residue = coeff of (x-a)^{-1}
    #   = (1/(k-1)!) d^{k-1}/dx^{k-1} [ num/dd ] |_{x=a}
    g_num = num
    g_den = dd
    # compute (k-1)-th derivative of num/dd at x=a via successive quotient-rule derivatives
    # represent h = g_num/g_den symbolically as a Qex element and differentiate.
    h = ctx.Qex(g_num) // ctx.Qex(g_den)
    for _ in 1:(k-1)
        h = derivative(h)
    end
    val = _eval_qex_at(ctx, h, aa)
    fact = factorial(big(k-1))
    return val // Qe(fact)
end

# evaluate a Qex element at x = aa (aa in Qe), returning Qe.
function _eval_qex_at(ctx::Ctx, f, aa)
    Qe = ctx.Qe
    num = numerator(f); den = denominator(f)
    nv = _eval_poly_at(num, aa, Qe)
    dv = _eval_poly_at(den, aa, Qe)
    return nv // dv
end

function _eval_poly_at(p, aa, Qe)
    s = Qe(0)
    d = degree(p)
    for i in 0:d
        s += coeff(p, i) * aa^i
    end
    return s
end

"Evaluate a Qex matrix at x = aa (Qe) -> Qe matrix."
function eval_x_at(ctx::Ctx, A, aa)
    m, n = size(A)
    Qe = ctx.Qe
    B = zero_matrix(Qe, m, n)
    a = Qe(aa)
    for i in 1:m, j in 1:n
        B[i,j] = _eval_qex_at(ctx, A[i,j], a)
    end
    return B
end

"""
    residue_at_infinity(ctx, A) -> matrix over Qe

Residue at x=∞ of A dx, computed DIRECTLY from the leading 1/x behaviour:
    R∞ = - lim_{x→∞} x·A(x) = -(coeff of 1/x in A) .
For a Fuchsian system this equals -(sum of ALL finite residues), but the previous
implementation summed only the *rational* poles (output of `singular_points_x`), silently
dropping irreducible-quadratic factors and giving a wrong R∞ when algebraic letters are
present (e.g. hexabox sec219's deg-2 letters).  Computing the limit directly is exact and
agnostic to the denominator factorization.
"""
function residue_at_infinity(ctx::Ctx, A)
    m, n = size(A)
    Qe = ctx.Qe
    S = zero_matrix(Qe, m, n)
    for i in 1:m, j in 1:n
        f = A[i,j]; iszero(f) && continue
        num = numerator(f); den = denominator(f)
        dn = degree(num); dd = degree(den)
        # coeff of 1/x in num/den: nonzero only when dd - dn == 1 (entry is O(1/x) and ≠ o(1/x))
        if dd - dn == 1
            S[i,j] = - leading_coefficient(num) // leading_coefficient(den)
        elseif dd - dn <= 0
            # not Fuchsian at ∞ — fall back to the (still well-defined) Laurent route via x→1/w
            g = subst_x_inv(ctx, f)            # f(1/w)
            h = - g // ctx.x^2                 # -f(1/w)/w² ; residue at w=0
            S[i,j] = _residue_scalar(ctx, h, Qe(0))
        end
    end
    return S
end

"""
    poly_pole_order(ctx, A) -> max pole order in x over all entries (1 = Fuchsian-candidate)

Returns the maximum multiplicity of any (finite) denominator factor, AND a flag for pole at
infinity higher than simple.  Used by is_fuchsian.
"""
function poly_pole_order(ctx::Ctx, A)
    maxord = 0
    m, n = size(A)
    for i in 1:m, j in 1:n
        a = A[i,j]
        iszero(a) && continue
        d = denominator(a)
        is_constant(d) && continue
        fd = factor(d)
        for (_, e) in fd
            maxord = max(maxord, e)
        end
    end
    return maxord
end

"""
    pole_profile(ctx, A) -> Dict{Any,Int}

Per-singularity pole order of the connection matrix A in x.  Keys are the finite singular
points (as Qe elements) and the symbol `:inf`.  Order = max over entries of the (x-a)
denominator multiplicity (finite) resp. 2 + max(deg num - deg den) (infinity; Fuchsian ⟺ 1).
A Fuchsian system has every value ≤ 1.
"""
function pole_profile(ctx::Ctx, A)
    m, n = size(A)
    prof = Dict{Any,Int}()
    # finite poles
    for i in 1:m, j in 1:n
        a = A[i,j]; iszero(a) && continue
        d = denominator(a)
        is_constant(d) && continue
        for (p, e) in factor(d)
            if degree(p) == 1
                c1 = coeff(p,1); c0 = coeff(p,0)
                root = -c0 // c1
                # true multiplicity in the reduced fraction (Nemo keeps num/den coprime,
                # so e is already the pole order of this entry at `root`)
                prof[root] = max(get(prof, root, 0), e)
            else
                # algebraic factor — record under the polynomial itself
                prof[p] = max(get(prof, p, 0), e)
            end
        end
    end
    # infinity: order = 2 + max_{ij}(deg num - deg den); ≤ 1 ⟺ A = O(1/x)
    dinf = -10^9
    for i in 1:m, j in 1:n
        a = A[i,j]; iszero(a) && continue
        dinf = max(dinf, degree(numerator(a)) - degree(denominator(a)))
    end
    prof[:inf] = dinf == -10^9 ? 0 : dinf + 2
    return prof
end

"""
    laurent_coeff_matrix(ctx, A, a, k) -> Qe-matrix

Coefficient of (x-a)^{-k} in the Laurent expansion of A at the finite point x=a (a ∈ Qe),
for k ≥ 1.  Computed exactly via  (1/(ord-k)!) d^{ord-k}/dx^{ord-k} [(x-a)^ord · A_{ij}] |_{x=a}
where ord is the entry's pole order.  Entries regular at a contribute 0.
"""
function laurent_coeff_matrix(ctx::Ctx, A, a, k::Int)
    m, n = size(A)
    Qe = ctx.Qe
    aa = Qe(a)
    L = zero_matrix(Qe, m, n)
    x = ctx.x
    Px = parent(numerator(x))
    lin = gen(Px) - Px(aa)
    for i in 1:m, j in 1:n
        f = A[i,j]; iszero(f) && continue
        # pole multiplicity μ of this entry at a
        dd = denominator(f); μ = 0
        while true
            q, r = divrem(dd, lin)
            iszero(r) || break
            dd = q; μ += 1
        end
        μ < k && continue
        # coeff of (x-a)^{-k}: (1/(μ-k)!) d^{μ-k} [(x-a)^μ f] at x=a
        h = (x - ctx.Qex(aa))^μ * f
        for _ in 1:(μ-k)
            h = derivative(h)
        end
        L[i,j] = _eval_qex_at(ctx, h, aa) // Qe(factorial(big(μ-k)))
    end
    return L
end

"""
    subst_x_inv(ctx, f) -> Qex element  f(1/x)

Exact substitution x → 1/x on a scalar rational function in Qex.
"""
function subst_x_inv(ctx::Ctx, f)
    iszero(f) && return ctx.Qex(0)
    num = numerator(f); den = denominator(f)
    Px = parent(num)
    dn = degree(num); dd = degree(den)
    # f(1/x) = (Σ n_i x^{-i}) / (Σ d_i x^{-i}) = (Σ n_i x^{D-i}) / (Σ d_i x^{D-i})  with D = max(dn,dd)
    D = max(dn, dd)
    revn = Px([ (D-i) <= dn ? coeff(num, D-i) : zero(ctx.Qe) for i in 0:D ])
    revd = Px([ (D-i) <= dd ? coeff(den, D-i) : zero(ctx.Qe) for i in 0:D ])
    return ctx.Qex(revn) // ctx.Qex(revd)
end

"Entrywise x → 1/x on a Qex matrix."
function subst_x_inv_mat(ctx::Ctx, A)
    m, n = size(A)
    B = zero_matrix(ctx.Qex, m, n)
    for i in 1:m, j in 1:n
        B[i,j] = subst_x_inv(ctx, A[i,j])
    end
    return B
end

"""
    is_fuchsian(ctx, A) -> Bool

A is Fuchsian in x iff every finite singularity is a simple pole AND the pole at infinity is
simple.  Pole at infinity simple <=> A = O(1/x) as x->∞ entrywise (deg den >= deg num + 1).
"""
function is_fuchsian(ctx::Ctx, A)
    # finite poles simple?
    if poly_pole_order(ctx, A) > 1
        return false
    end
    # behaviour at infinity: each entry must be O(1/x)
    m, n = size(A)
    for i in 1:m, j in 1:n
        a = A[i,j]
        iszero(a) && continue
        if degree(numerator(a)) >= degree(denominator(a))
            return false
        end
    end
    return true
end

end # module

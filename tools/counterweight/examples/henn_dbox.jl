# henn_dbox.jl — the massless planar double box of Henn arXiv:1304.1806 (eq. derivativeint).
#
# Canonical (eps-form) connection in the ratio variable x = t/s:
#     ∂_x f = eps [ a/x + b/(1+x) ] f
# with the constant 8×8 rational matrices a, b given verbatim in 1304.1806 (transcribed below
# from refs/1304.1806/qcdintegrals_arxiv2.tex).  Singular points {0, -1, ∞}; alphabet {x, 1+x}.
#
# This module exposes:
#   henn_a(), henn_b()         — the constant QQ matrices.
#   henn_canonical(ctx)        — the eps-form connection A = eps[a/x + b/(1+x)] as a Qex matrix.
#   henn_scramble(ctx, A; ...) — apply a nontrivial eps&x-dependent gauge to de-canonicalize A
#                                into a generic (non-UT) basis, simulating a raw Kira basis.

module HennDbox

using Nemo
using ..RatFunc

export henn_a, henn_b, henn_canonical, henn_scramble

function henn_a()
    # rows of matrix a (Henn 1304.1806)
    rows = [
        [-2, 0, 0, 0, 0, 0, 0, 0],
        [ 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0, 0, 0, 0, 0, 0, 0, 0],
        [ QQ(3,2), 0, 0, 0, -2, 0, 0, 0],
        [ QQ(-1,2), QQ(1,2), 0, 0, 0, -2, 0, 0],
        [ -3, -3, 0, 0, 4, 12, -2, 0],
        [ QQ(9,2), 3, -3, -1, -4, -18, 1, 1],
    ]
    M = zero_matrix(QQ, 8, 8)
    for i in 1:8, j in 1:8
        M[i,j] = QQ(rows[i][j])
    end
    return M
end

function henn_b()
    rows = [
        [ 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0, 0, 0, 0, 0, 0, 0, 0],
        [ 0, 0, 0, 0, 0, 0, 0, 0],
        [ QQ(-3,2), 0, 3, 0, 1, 0, 0, 0],
        [ 0, 0, 0, 0, 0, 2, 0, 0],
        [ 3, 6, 6, 2, -4, -12, 2, 2],
        [ QQ(-9,2), -3, 3, -1, 4, 18, -1, -1],
    ]
    M = zero_matrix(QQ, 8, 8)
    for i in 1:8, j in 1:8
        M[i,j] = QQ(rows[i][j])
    end
    return M
end

"Canonical eps-form connection A(eps,x) = eps[a/x + b/(1+x)] as an 8×8 Qex matrix."
function henn_canonical(ctx::Ctx)
    a = henn_a(); b = henn_b()
    e = ctx.eps; x = ctx.x
    A = zero_matrix(ctx.Qex, 8, 8)
    for i in 1:8, j in 1:8
        A[i,j] = e*( ctx.Qex(a[i,j])//x + ctx.Qex(b[i,j])//(x+1) )
    end
    return A
end

"""
    henn_scramble(ctx, A; mode=:lowertri)

De-canonicalize A by a nontrivial gauge B(eps,x), producing A' = B A B^{-1} + (dB)B^{-1} that is
Fuchsian but NOT eps-form — a stand-in for a generic raw IBP basis.  The gauge mixes UT masters
with eps-dependent, x-dependent rational coefficients (eps^0 eigenvalue shifts + polynomial
off-diagonal mixing), exactly the kind of non-UT contamination the engine must undo.
"""
function henn_scramble(ctx::Ctx, A; verbose=false)
    e = ctx.eps; x = ctx.x
    n = size(A,1)
    I = identity_matrix(ctx.Qex, n)
    # Build a composite gauge:
    #  (1) integer eps^0 eigenvalue shift via a diagonal balance B1 on master 6 between poles 0 and -1
    #  (2) a constant invertible mixing C0 (unimodular integer matrix) that scrambles the basis
    #  (3) an x-dependent shear adding (poly in x) * (lower master) into an upper master
    # All keep the system Fuchsian with the same singular set {0,-1,inf}.
    # An x-INDEPENDENT (constant in x) gauge C(eps) preserves Fuchsianity exactly (dC/dx = 0),
    # and when C depends on eps it breaks eps-form — this is precisely how a raw Kira/IBP master
    # basis relates to the UT basis: raw masters are eps-dependent rational combinations of UT
    # masters.  We build a nontrivial such C(eps):
    #   (1) an eps-dependent diagonal rescaling D(eps) (shifts the ε-weight grading away from UT)
    #   (2) a constant unimodular integer mixing C0 (basis rotation among masters)
    #   (3) an eps-dependent lower-triangular shear (adds (1+eps)*lower master into an upper one)
    D = identity_matrix(ctx.Qex, n)
    D[1,1] = ctx.Qex(1+e)
    D[6,6] = ctx.Qex(1)//(1+2*e)
    A1 = _gauge(ctx, A, D)
    C0 = identity_matrix(ctx.Qex, n)
    C0[5,1] = ctx.Qex(2); C0[7,3] = ctx.Qex(-1); C0[8,2] = ctx.Qex(3); C0[6,5] = ctx.Qex(1)
    C0[4,1] = ctx.Qex(-2)
    A2 = _gauge(ctx, A1, C0)
    Sh = identity_matrix(ctx.Qex, n)
    Sh[7,6] = ctx.Qex(1+e)
    Sh[8,1] = ctx.Qex(e)
    A3 = _gauge(ctx, A2, Sh)
    return A3, (Sh*C0*D)
end

function _gauge(ctx::Ctx, A, T)
    Tinv = inv(T)
    dT = zero_matrix(ctx.Qex, size(T)...)
    m,n = size(T)
    for i in 1:m, j in 1:n
        dT[i,j] = derivative(T[i,j])
    end
    return T*A*Tinv + dT*Tinv
end

end # module

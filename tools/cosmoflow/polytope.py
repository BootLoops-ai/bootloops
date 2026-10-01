#!/usr/bin/env python3
r"""
polytope.py -- Cosmological one-loop n-site (n-gon) Baikov polynomial and
linear polytope-facet forms q_G.

Geometry: for an n_s-site 1-loop FRW graph, the loop-tetrahedron generalises
to an (n_s)-simplex with base vertices V_1..V_{n_s} (external kinematics) and
loop apex L.  Squared edge lengths:
    base:   d(V_i,V_j)^2 = X[i,j]^2    (n_s*(n_s-1)/2 externals)
    apex:   d(L,V_i)^2   = y[i]^2      (n_s loop lengths, y_e >= 0)
The Baikov polynomial is the (n_s+2)x(n_s+2) Cayley--Menger determinant on
these squared lengths (arXiv:2408.16386 eq. mCM / eq. ukchi):
    B(y;X) = CM_{n_s+1}(d^2)  =  (-1)^{n_s} 2^{n_s} (n_s!)^2 Vol^2.
In d = n_s spatial dimensions, gamma = (d - n_s - 1)/2 = -1/2 so the twist
is mu_d ~ B^{-1/2}.  B is quadratic in each z_e = y_e^2.

The linear-denominator "cosmological polytope" facet forms q_G are linear in
the y_e (NOT the z_e); for the sub-sector G_{j,j+1} of the n-gon
    q_{G_{j,j+1}} = sum_s X_s + c * y_{j,j+1},   c = 2 (DERIVED),
c=2 is DERIVED, not assumed -- three independent routes
(from-scratch flat-space wavefunction pole census, canonical-form ratio test,
vertex-facet geometry from the 1709.02813 vertex definition) agree exactly;
the c=1 printed in 2408.16386's triangle section is a typo internally
contradicted three times in the same paper.  c=1 stays reachable by
explicit argument (L2 and the 16x16 connection are proven c-independent).

n_s = 3 (tetrahedron) is the tested/validated special case.
n_s = 2, >=4 are constructed by the same CM formula but NOT independently
validated -- use for structure, not for gated numerics.
"""
import sympy as sp


def cayley_menger(d2):
    """Cayley--Menger determinant on an (n x n) squared-distance matrix d2
    (sympy Matrix, d2[i,i]=0, d2 symmetric).  Returns det of the bordered
    (n+1)x(n+1) CM matrix."""
    n = d2.shape[0]
    M = sp.zeros(n + 1, n + 1)
    for i in range(n):
        M[0, i + 1] = M[i + 1, 0] = 1
        for j in range(n):
            M[i + 1, j + 1] = d2[i, j]
    return M.det()


def baikov_B(n_s, y_list, X_list):
    """Baikov polynomial B = CM_{n_s+1} of the loop simplex.
        y_list : length-n_s iterable of apex-to-base edge lengths y_i
        X_list : length-n_s*(n_s-1)//2 iterable of base edge lengths X_{ij},
                 upper-triangular row-major order (i<j): X_{12},X_{13},...,X_{1n},
                 X_{23},...,X_{n-1,n}.
    Returns a sympy expression, quadratic in each y_i^2."""
    y_list = list(y_list); X_list = list(X_list)
    if len(y_list) != n_s:
        raise ValueError(f"need {n_s} y-lengths, got {len(y_list)}")
    if len(X_list) != n_s * (n_s - 1) // 2:
        raise ValueError(f"need {n_s*(n_s-1)//2} X-lengths, got {len(X_list)}")
    n = n_s + 1                       # points: V_1..V_{n_s}, L
    d2 = sp.zeros(n, n)
    k = 0
    for i in range(n_s):
        for j in range(i + 1, n_s):
            d2[i, j] = d2[j, i] = X_list[k] ** 2
            k += 1
    for i in range(n_s):
        d2[i, n_s] = d2[n_s, i] = y_list[i] ** 2
    return sp.expand(cayley_menger(d2))


def q_subgraph(subgraph_spec, y_list, X_list, c=2):
    """Linear polytope-facet form q_G.
        subgraph_spec : (y_indices, X_indices, c)  OR  a single int j
                        (shorthand for the G_{j,j+1} facet: all X's + c*y[j]).
        c             : coefficient used by the int shorthand ONLY (a tuple
                        spec carries its own c).  Default c=2 -- the deleted
                        loop edge is cut TWICE (derived; see module
                        docstring).  Pass c=1 explicitly for the superseded
                        2408.16386 print convention.
    Returns  sum(X_list[i] for i in X_indices) + c * sum(y_list[i] for i in y_indices)."""
    if isinstance(subgraph_spec, int):
        y_idx, X_idx = (subgraph_spec,), tuple(range(len(X_list)))
    else:
        y_idx, X_idx, c = subgraph_spec
    return (sum(X_list[i] for i in X_idx)
            + c * sum(y_list[i] for i in y_idx))


# --- n_s = 3 tetrahedron: the validated special case, used by the WORKED
# --- EXAMPLE (examples/triangle: oracle, maxcut) ---------------------------
# Edge convention:
#   base vertices A,B,C = V_1,V_2,V_3 ;  apex L
#   |AB|=X2, |AC|=X1, |BC|=X3    (=> X_list order (X_{12},X_{13},X_{23})=(X2,X1,X3))
#   |LA|=y12, |LB|=y23, |LC|=y31 (=> y_list order (y_1,y_2,y_3)=(y12,y23,y31))
# so opposite pairs are (y12,X3),(y23,X1),(y31,X2).
y12, y23, y31, X1, X2, X3 = sp.symbols('y12 y23 y31 X1 X2 X3')
z1, z2, z3 = sp.symbols('z1 z2 z3')  # z_e = y_e^2


def baikov_B3(y12, y23, y31, X1, X2, X3):
    """n_s=3 tetrahedron Baikov polynomial (= 288*Vol^2).  Validated."""
    return baikov_B(3, (y12, y23, y31), (X2, X1, X3))


def q_G12(y12, X1, X2, X3, c=2):
    """q_{G_{12}} = X1+X2+X3 + c*y12.

    c=2 is the DERIVED physical coefficient (three independent exact
    routes); the c=1 printed in 2408.16386's triangle section is a typo,
    internally contradicted by its own bubble section, general marking
    rule, and figure caption.  c=1 remains reachable by explicit argument
    (L2/periods proven c-independent)."""
    return X1 + X2 + X3 + c * y12


# Precomputed expansion in Baikov vars z_e = y_e^2 (polynomial, deg<=2 each).
B_z = sp.expand(baikov_B(3, (sp.sqrt(z1), sp.sqrt(z2), sp.sqrt(z3)),
                         (X2, X1, X3)))
B_num = sp.lambdify((z1, z2, z3, X1, X2, X3), B_z, modules='mpmath')


if __name__ == '__main__':
    # Regular unit tetrahedron -> Vol^2 = 1/72 -> CM = 4.
    print("CM(unit tetra) =", baikov_B3(1, 1, 1, 1, 1, 1), "  (expect 4)")
    for zv in (z1, z2, z3):
        print(f"  deg_{zv}(B) =", sp.degree(B_z, zv), "  (expect 2)")
    # n_s=2 (triangle simplex): regular unit -> CM_3 = -3 (2*Area^2*16 sign conv).
    print("CM(n_s=2, unit) =", baikov_B(2, (1, 1), (1,)))
    # n_s=4 sanity: quadratic in each y.
    ys = sp.symbols('y1:5'); Xs = sp.symbols('X1:7')
    B4 = baikov_B(4, ys, Xs)
    print("  deg_{y1^2}(B4) =", sp.degree(sp.Poly(B4, ys[0]**2), ys[0]**2))

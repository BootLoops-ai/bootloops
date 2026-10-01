# SPDX-License-Identifier: GPL-2.0-or-later
#
# Derived from SageMath (https://www.sagemath.org),
# sage/schemes/hyperelliptic_curves/mestre.py (Florian Bouyer, Marco Streng;
# Copyright (C) 2011, 2012, 2013 Florian Bouyer, Marco Streng; 2025 Sabrina
# Kunzweiler, Gareth Ma, Giacomo Pope) and
# sage/schemes/hyperelliptic_curves/invariants.py (Nick Alexander; Copyright (C) 2008
# Nick Alexander; 2025 Sabrina Kunzweiler, Gareth Ma, Giacomo Pope), licensed under the
# GNU General Public License v2.0 or later. This sympy transliteration of Mestre's
# conic-and-cubic construction and of the Clebsch and Igusa-Clebsch invariants
# (J.-F. Mestre, "Construction de courbes de genre 2 a partir de leurs modules", 1991;
# K. Lauter and T. Yang, J. Number Theory 131 (2011) 936; invariants after J. Igusa,
# 1960) keeps SageMath's conventions and normalizations so that its published doctest
# values serve as exact self-checks; the transvectant routine is written independently.
# It is Copyright (c) 2026 Anthropic, PBC; created by Matthew D. Schwartz, code written
# by Claude (Anthropic) under his supervision, and is distributed under the same GNU
# General Public License v2.0 or later. The rest of this package is MIT-licensed; see
# ../NOTICE (tools/eichler/NOTICE).
#
# mestre_port.py: Mestre/Igusa-Clebsch toolkit for eichler.genus2 (no Sage dependency).
import sympy as sp
from functools import lru_cache

X, Y = sp.symbols('mx my')

def transvect(f, g, k, m=None, n=None):
    """(f g)_k for homogeneous binary forms f (deg m), g (deg n) in (X, Y)."""
    if m is None: m = sp.total_degree(sp.Poly(f, X, Y))
    if n is None: n = sp.total_degree(sp.Poly(g, X, Y))
    pref = sp.Rational(sp.factorial(m-k)*sp.factorial(n-k), sp.factorial(m)*sp.factorial(n))
    tot = sp.Integer(0)
    for j in range(k+1):
        tot += (-1)**j * sp.binomial(k, j) * sp.diff(f, X, k-j, Y, j) * sp.diff(g, X, j, Y, k-j)
    return sp.expand(pref*tot)

def clebsch_ABCD(f6):
    """Clebsch A,B,C,D of a homogeneous binary sextic f6(X,Y)."""
    f = sp.expand(f6)
    i_ = transvect(f, f, 4, 6, 6)                    # deg 4
    Delta = transvect(i_, i_, 2, 4, 4)               # deg 4
    y1 = transvect(f, i_, 4, 6, 4)                   # deg 2
    y2 = transvect(i_, y1, 2, 4, 2)                  # deg 2
    y3 = transvect(i_, y2, 2, 4, 2)                  # deg 2
    A = transvect(f, f, 6, 6, 6)
    B = transvect(i_, i_, 4, 4, 4)
    C = transvect(i_, Delta, 4, 4, 4)
    D = transvect(y3, y1, 2, 2, 2)
    return tuple(sp.expand(t) for t in (A, B, C, D))

def clebsch_to_igusa(A, B, C, D):
    I2 = -120*A
    I4 = -720*A**2 + 6750*B
    I6 = 8640*A**3 - 108000*A*B + 202500*C
    I10 = (-62208*A**5 + 972000*A**3*B + 1620000*A**2*C
           - 3037500*A*B**2 - 6075000*B*C - 4556250*D)
    return (I2, I4, I6, I10)

def igusa_clebsch(f6):
    return clebsch_to_igusa(*clebsch_ABCD(f6))

def mestre_xyz(I2, I4, I6, I10):
    x = 8*(1 + 20*I4/(I2**2))/225
    y = 16*(1 + 80*I4/(I2**2) - 600*I6/(I2**3))/3375
    z = -64*(-10800000*I10/(I2**5) - 9 - 700*I4/(I2**2) + 3600*I6/(I2**3)
             + 12400*I4**2/(I2**4) - 48000*I4*I6/(I2**5))/253125
    return sp.together(x), sp.together(y), sp.together(z)

def mestre_conic_matrix(x, y, z):
    return sp.Matrix([
        [x + 6*y, 6*x**2 + 2*y, 2*z],
        [6*x**2 + 2*y, 2*z, 9*x**3 + 4*x*y + 6*y**2],
        [2*z, 9*x**3 + 4*x*y + 6*y**2, 6*x**2*y + 2*y**2 + 3*x*z]])

def mestre_cijk(x, y, z):
    c = {}
    c[(1,1,1)] = 12*x*y - 2*y/3 - 4*z
    c[(1,1,2)] = -18*x**3 - 12*x*y - 36*y**2 - 2*z
    c[(1,1,3)] = -9*x**3 - 36*x**2*y - 4*x*y - 6*x*z - 18*y**2
    c[(1,2,2)] = c[(1,1,3)]
    c[(1,2,3)] = -54*x**4 - 36*x**2*y - 36*x*y**2 - 6*x*z - 4*y**2 - 24*y*z
    c[(1,3,3)] = (-sp.Rational(27,2)*x**4 - 72*x**3*y - 6*x**2*y - 9*x**2*z
                  - 39*x*y**2 - 36*y**3 - 2*y*z)
    c[(2,2,2)] = -27*x**4 - 18*x**2*y - 6*x*y**2 - sp.Rational(8,3)*y**2 + 2*y*z
    c[(2,2,3)] = 9*x**3*y - 27*x**2*z + 6*x*y**2 + 18*y**3 - 8*y*z
    c[(2,3,3)] = (-sp.Rational(81,2)*x**5 - 27*x**3*y - 9*x**2*y**2 - 4*x*y**2
                  + 3*x*y*z - 6*z**2)
    c[(3,3,3)] = (sp.Rational(27,2)*x**4*y - sp.Rational(27,2)*x**3*z + 9*x**2*y**2
                  + 3*x*y**3 - 6*x*y*z + sp.Rational(4,3)*y**3 - 10*y**2*z)
    return c

def curve_from_parametrization(x, y, z, F1, F2, F3, tvar):
    c = mestre_cijk(x, y, z)
    F = {1: F1, 2: F2, 3: F3}
    f = sp.Integer(0)
    for (i, j, k), cv in c.items():
        # symmetric sum convention: Sage sums each listed (i<=j<=k) ONCE with its c
        f += cv * F[i]*F[j]*F[k]
    return sp.expand(f)

if __name__ == '__main__':
    # verification vs Sage doctest vectors
    t = sp.Symbol('t')
    f1 = X**6 + Y**6
    print('clebsch(x^6+y^6):', clebsch_ABCD(f1), ' expect (2, 2/3, -2/9, 0)')
    print('IC(x^6+y^6):', igusa_clebsch(f1), ' expect (-240, 1620, -119880, -46656)')
    p = t**6 + t**5 + t**4 + t**2 + 2
    ph = sp.expand(sum(sp.Poly(p, t).coeff_monomial(t**i)*X**i*Y**(6-i) for i in range(7)))
    print('IC(x^6+x^5+x^4+x^2+2):', igusa_clebsch(ph), ' expect (-496, 6220, -955932, -1111784)')
    # Mestre conic check vs Sage: Mestre_conic([1,2,3,4]) coefficients
    x_, y_, z_ = mestre_xyz(1, 2, 3, 4)
    L = mestre_conic_matrix(x_, y_, z_)
    den = sp.lcm([sp.denom(sp.together(L[i, j])) for i in range(3) for j in range(3)])
    Lc = sp.expand(L*den)
    # conic form: u^T L u: coefficient of u^2 = L00, uv = 2*L01 ...
    print('conic check [1,2,3,4]: u2:', Lc[0,0], ' uv:', 2*Lc[0,1], ' v2:', Lc[1,1],
          ' uw:', 2*Lc[0,2], ' vw:', 2*Lc[1,2], ' w2:', Lc[2,2])
    print(' expect ratios of (-2572155000, -317736000, 1250755459200, 2501510918400, 39276887040, 2736219686912)')

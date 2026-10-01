# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
r"""
feynman.py -- build the arbitrary-precision integrand for a finite (eps^0)
scalar Feynman integral from its propagator list, with automatic closed-form
elimination of every Feynman parameter that enters F linearly.

Pipeline
--------
1.  Symanzik U, F via pySecDec's LoopIntegralFromPropagators (sympy polynomials).
2.  Assert the integral is finite at eps^0 (the caller is responsible; this tool
    targets the eps^0 coefficient I = Gamma(N - L d/2) * \int_simplex U^a / F^b,
    with a = N-(L+1)d/2, b = N-Ld/2 at eps=0).
3.  Cheng-Wu: set the last parameter x_{n-1}=1, integrate the remaining over
    [0,inf)^{n-1}.
4.  ANALYTIC REDUCTION: iteratively integrate out variables linear in the
    denominator (a==0, the finite-box case), through two closed-form stages:
      stage 1  \int_0^inf dv / (A v + B)^2         = 1/(A B)
      stage 2  \int_0^inf dv / ((A1 v+A0)(B1 v+B0))
               = log(A1 B0 / (A0 B1)) / (A1 B0 - A0 B1)
    (all coefficient polynomials positive in the Euclidean region).  Each
    stage removes one numerical dimension exactly; `reduce_linear_full` runs
    the iteration, `reduce_linear` exposes stage 1's factor form.
5.  Emit a gmpy2-lambdified integrand factory over the surviving [0,inf)^dim.

The closed-form reduction targets the finite a==0, b==2 case
(N - (L+1)d/2 = 0 and N - L d/2 = 2 at eps=0, d=4), which covers the
one-loop box, the 3-loop self-energy benchmark, and the generic "two more
propagators than 2L" finite topology.  For other (a,b) the reduction step is
disabled and the full (n-1)-dim integrand U^a/F^b is emitted (finiteness of
the target order stays the caller's responsibility).
"""
import sympy as sp


PROP_CACHE = {}


def build_UF(propagators, loop_momenta, external_momenta, replacement_rules):
    """Return (U, F, xs, params) as sympy objects.  params = (s,t,m2,...) symbols
    found in the replacement rules / propagators."""
    from pySecDec.loop_integral import LoopIntegralFromPropagators
    li = LoopIntegralFromPropagators(
        propagators=propagators,
        loop_momenta=loop_momenta,
        external_momenta=external_momenta,
        replacement_rules=replacement_rules,
        regulators=['eps'],
    )
    U = sp.sympify(str(li.U))
    F = sp.sympify(str(li.F))
    xs = sp.symbols(' '.join(str(x) for x in li.Feynman_parameters))
    if not isinstance(xs, (list, tuple)):
        xs = (xs,)
    L = int(li.L)
    Nprop = len(propagators)
    return U, F, list(xs), L, Nprop


def reduce_linear(F, xs, b=2):
    """Cheng-Wu (x_last = 1), then the stage-1 closed-form elimination: if some
    variable v is linear in the (single, squared) denominator polynomial,
    integrate it out via \\int_0^inf dv/(A v + B)^2 = 1/(A B).

    Returns (den_factors, surviving_vars): den_factors is [Fc] (denominator
    Fc^2, no elimination possible) or [A, B] (denominator A*B, each factor to
    power 1, one variable eliminated).  After stage 1 the denominator is a
    product of two power-1 factors, so a further variable can be eliminated
    only by partial fractions -- that continuation is `reduce_linear_full`,
    which callers wanting the maximal reduction should use directly.
    """
    if b != 2:
        raise NotImplementedError("reduce_linear implemented for denominator power b=2")
    x = list(xs)
    Fc = sp.expand(F.subs(x[-1], 1))
    remaining = x[:-1]
    den_factors = [Fc]  # denominator = prod(den_factors)^2 while len == 1
    for v in list(remaining):
        if len(den_factors) == 1:
            Fd = sp.expand(den_factors[0])
            if sp.degree(Fd, v) == 1:
                A = Fd.coeff(v, 1)
                B = Fd.coeff(v, 0)
                # \int_0^inf dv/(A v + B)^2 = 1/(A B)
                den_factors = [sp.expand(A), sp.expand(B)]
                remaining.remove(v)
    return den_factors, remaining


def reduce_linear_full(F, xs, b=2):
    """Cheng-Wu (x_last = 1), then iterated closed-form elimination of
    variables linear in the denominator, through the two implemented stages
    (up to two exact dimension drops for b=2).

    Stage 1 (via `reduce_linear`): a variable linear in the squared
    denominator polynomial:  \\int_0^inf dv/(A v + B)^2 = 1/(A B).
    Stage 2: a variable linear in BOTH surviving power-1 factors -- partial
    fractions give
        \\int_0^inf dv/((A1 v + A0)(B1 v + B0))
            = log(A1 B0 / (A0 B1)) / (A1 B0 - A0 B1),
    with the removable point A1 B0 == A0 B1 (proportional factors) evaluating
    to 1/(A0 B1).  Validity as everywhere in this module: every coefficient
    polynomial positive on the open domain (deep-Euclidean region).

    Each stage removes one numerical dimension exactly.  A variable linear in
    only ONE power-1 factor is never eliminated (that integral diverges), and
    a variable of higher degree in either factor is left for
    `reduce_one_quadratic` or numerics.  After stage 2 the integrand is
    transcendental (a log ratio), which ends the iteration.

    Returns (integrand_expr, surviving_vars): the exact integrand over
    [0,inf)^len(surviving_vars) with all closed-form factors included.
    """
    if b != 2:
        raise NotImplementedError("reduce_linear_full implemented for denominator power b=2")
    den_factors, remaining = reduce_linear(F, xs, b=b)
    if len(den_factors) == 1:
        return 1/den_factors[0]**2, remaining
    A, B = den_factors
    for v in list(remaining):
        if sp.degree(A, v) == 1 and sp.degree(B, v) == 1:
            A1, A0 = A.coeff(v, 1), A.coeff(v, 0)
            B1, B0 = B.coeff(v, 1), B.coeff(v, 0)
            num = sp.expand(A1*B0)
            den = sp.expand(B1*A0)
            remaining.remove(v)
            if sp.expand(num - den) == 0:
                expr = 1/(A0*B1)
            else:
                expr = sp.Piecewise((1/(A0*B1), sp.Eq(num, den)),
                                    (sp.log(num/den)/(num - den), True))
            return expr, remaining
    return 1/(A*B), remaining


def reduce_one_quadratic(F, xs, var, b=2):
    """Integrate ONE Feynman parameter `var` that appears quadratically in F out of
    1/F^b over [0,inf), in closed form (the "box" reduction).  Returns the new
    integrand expression (a sympy expr in the remaining variables, generally
    containing atan/atanh/sqrt) and the list of remaining variables.

    For b=2:  \\int_0^inf dx/(P x^2 + Q x + R)^2
              = -Q/(D R) + (2P/D) J1 ,   D = 4PR - Q^2 ,
              J1 = (2/sqrt(D))(pi/2 - atan(Q/sqrt(D)))         if D>0
                 = -(1/sqrt(-D)) log((Q-sqrt(-D))/(Q+sqrt(-D)))  if D<0
    valid for P>0 and the quadratic positive on [0,inf) (Euclidean region).

    This is the symbolic generalisation of bench_box._I2 -- use it to collapse one
    loop parameter of a box-like sub-integral before handing the rest to tanh-sinh.
    The branch (D>0 vs D<0) is resolved numerically at evaluation time, so the
    returned expression keeps both via a sympy Piecewise.
    """
    if b != 2:
        raise NotImplementedError("reduce_one_quadratic implemented for b=2")
    Fc = sp.expand(F.subs(xs[-1], 1)) if xs[-1] in F.free_symbols else sp.expand(F)
    P = Fc.coeff(var, 2)
    Q = Fc.coeff(var, 1)
    R = Fc.coeff(var, 0)
    if P == 0:
        raise ValueError("variable %s is not quadratic in F" % var)
    D = sp.expand(4*P*R - Q**2)
    sD = sp.sqrt(sp.Abs(D))
    J1_pos = (2/sD)*(sp.pi/2 - sp.atan(Q/sD))
    J1_neg = -(1/sD)*sp.log((Q - sD)/(Q + sD))
    J1 = sp.Piecewise((J1_pos, D > 0), (J1_neg, True))
    expr = (-Q/(D*R)) + (2*P/D)*J1
    remaining = [v for v in xs[:-1] if v != var]
    return expr, remaining


def assemble_integrand_expr(U, F, xs, L, Nprop, param_values, d=4):
    """Assemble the eps^0 integrand from already-built Symanzik polynomials.

    Returns (expr, surviving_vars) for the integrand over
    [0,inf)^len(surviving_vars), with numeric param_values (dict
    sympy-symbol -> sympy number) substituted.

    a = N-(L+1)d/2 == 0 and b = N-Ld/2 == 2: the maximal closed-form linear
    reduction runs (`reduce_linear_full` -- up to two exact dimension drops).
    Any other (a, b): the reduction step is disabled and the full
    (n-1)-dim Cheng-Wu integrand U^a/F^b is emitted unreduced; finiteness of
    the target order is the caller's responsibility (this mirrors the
    projective pathway of parametric.py, which evaluates such integrands
    directly).
    """
    a = Nprop - (L+1)*d//2   # exponent of U at eps=0
    b = Nprop - L*d//2       # exponent of F at eps=0
    if a == 0 and b == 2:
        Fnum = F.subs(param_values)
        expr, remaining = reduce_linear_full(Fnum, xs, b=2)
        return sp.expand(expr), remaining
    Ucw = sp.expand(U.subs(param_values).subs(xs[-1], 1))
    Fcw = sp.expand(F.subs(param_values).subs(xs[-1], 1))
    expr = Ucw**a / Fcw**b
    return expr, list(xs[:-1])


def build_integrand_expr(propagators, loop_momenta, external_momenta,
                         replacement_rules, param_values):
    """Full convenience builder: returns (expr, surviving_vars) for the eps^0
    integrand over [0,inf)^len(surviving_vars), with numeric param_values
    (dict sympy-symbol -> sympy number) substituted.

    a=0, b=2 integrals get the maximal closed-form linear reduction; for any
    other (a, b) the full (n-1)-dim Cheng-Wu integrand U^a/F^b is emitted
    unreduced (see `assemble_integrand_expr`).
    """
    U, F, xs, L, Nprop = build_UF(propagators, loop_momenta,
                                  external_momenta, replacement_rules)
    return assemble_integrand_expr(U, F, xs, L, Nprop, param_values)


def make_factory(expr, remaining_vars):
    """Return an integrand_factory suitable for hiprec.integrate_qmc.

    The factory ignores its argument and returns a gmpy2-backed callable
    f(*x) over the `remaining_vars`.  Defined at module level via a closure
    builder so it survives multiprocessing fork (POSIX)."""
    import gmpy2
    src = str(expr)
    nvar = len(remaining_vars)
    varnames = [str(v) for v in remaining_vars]

    def factory(_):
        import sympy as _sp
        import gmpy2 as _g
        e = _sp.sympify(src)
        if nvar == 0:
            # fully reduced: the expression is an exact constant
            val = _g.mpfr(str(_sp.N(e, max(30, int(_g.get_context().precision/3.3)+5))))
            return lambda: val
        syms = _sp.symbols(' '.join(varnames))
        if nvar == 1:
            syms = (syms,)
        f = _sp.lambdify(syms, e,
                         modules=[{'sqrt': _g.sqrt, 'log': _g.log,
                                   'atan': _g.atan}, 'math'], cse=True)
        return f
    return factory, nvar

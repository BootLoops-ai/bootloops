"""mixalot.boxwalk.fiber — scalar-recurrence emission along the final ray.

The walk's final segment is a pure ray u_base + n*e_c: the fiber sequence
Z(u_base + n*e_c) mod p comes for free. Fit the dataset's own univariate
contiguity recurrence sum_j c_j(n) Z(n+j) = 0 with the annihilator engine
(find_recurrence_fast / nullspace_fast), 2-prime-stable (r,d), CRT-lift
coefficients, verify at a held-out prime.  Output: an exact integer scalar
operator per dataset — the refit trap broken BY the walk.

The fiber engine = dense oracle along the ray (always available; exact mod p).
Lifted operators carry a shift-algebra (Ore) calculus for reduction and
stitching: `op_mul`, `right_divides`, `gcrd`, `lclm`, `annihilates_mod`,
`gcrd_reduce` — fraction-free right pseudo-division, primitive GCRD by the
right Euclidean algorithm, LCLM by the extended algorithm, with exact
right-division certificates on every reduction.
"""
import importlib.util, math
import os as _os
from fractions import Fraction
import numpy as np
from . import core

# The recurrence fitter (find_recurrence_fast) is mixalot's own sha-pinned
# `annihilator` engine (mixalot/vendor/tools/annihilator.py, loaded through
# mixalot.engines.load -- source-hash-checked, fail-closed).  To run against
# another copy (e.g. the standalone tools/annihilator/annihilator.py), set
# BOXWALK_ANNIHILATOR=/path/to/annihilator.py.
_ANNI = _os.environ.get('BOXWALK_ANNIHILATOR') or ''


def _annihilator():
    if _ANNI:
        s = importlib.util.spec_from_file_location('annihilator', _ANNI)
        m = importlib.util.module_from_spec(s)
        s.loader.exec_module(m)
        return m
    from .. import engines
    return engines.load('annihilator')


def fiber_mod(spec, u_base, ray_class, depth, p):
    """Z(u_base + n*e_ray) mod p for n = 0..depth (dense-oracle engine)."""
    Q = core.dense_Q(spec, u_base, p)
    P = spec.polys[ray_class]
    Zs = []
    for _ in range(depth + 1):
        Zs.append(int(core.window_moments(Q, 0, p).ravel()[0]))
        Q = core.conv_mod(Q, P, p)
    return Zs


def emit_fiber_recurrence(spec, ray_class=None, u_base=None, depth=80,
                          rmax=20, dmax=12, primes=None, log=print):
    """Fit/lift/verify chain. Returns dict with (r, d), integer coeffs
    {(j,k): int}, and the held-out verification result."""
    anni = _annihilator()
    primes = primes or core.primes31(3)
    if ray_class is None:
        ray_class = max(spec.target_u, key=lambda k: spec.target_u[k])
    if u_base is None:
        u_base = {k: v for k, v in spec.target_u.items() if k != ray_class}
    fibers, hits = {}, {}
    for p in primes[:2]:
        fibers[p] = fiber_mod(spec, u_base, ray_class, depth, p)
        hits[p] = anni.find_recurrence_fast(fibers[p], p, rmax=rmax,
                                            dmax=dmax, prefer='ode_order')
        if hits[p] is None:
            return {'ok': False, 'error': f'no recurrence in range at p={p} '
                    f'(raise depth/rmax/dmax)'}
    (rA, dA, vA), (rB, dB, vB) = hits[primes[0]], hits[primes[1]]
    if (rA, dA) != (rB, dB):
        return {'ok': False, 'error': f'(r,d) unstable across primes: '
                f'{(rA, dA)} vs {(rB, dB)}'}
    r, d = rA, dA
    # normalize both nullvectors by the same pivot entry, CRT, ratrec
    j0 = next(i for i in range(len(vA)) if vA[i] % primes[0]
              and vB[i] % primes[1])
    pA, pB = primes[:2]
    nA = [v * pow(vA[j0], pA - 2, pA) % pA for v in vA]
    nB = [v * pow(vB[j0], pB - 2, pB) % pB for v in vB]
    M = pA * pB
    fracs = []
    for a, b in zip(nA, nB):
        rr, _ = core.crt_list([a, b], [pA, pB])
        fr = core.ratrec(rr, M)
        if fr is None:
            return {'ok': False, 'error': 'ratrec failed (coefficients too '
                    'large for 2 primes — extend lift)'}
        fracs.append(fr)
    lcm = 1
    for fr in fracs:
        lcm = lcm * fr.denominator // math.gcd(lcm, fr.denominator)
    ints = [int(fr * lcm) for fr in fracs]
    g = 0
    for v in ints:
        g = math.gcd(g, abs(v))
    ints = [v // max(g, 1) for v in ints]
    coeffs = {(j, k): ints[j * (d + 1) + k] for j in range(r + 1)
              for k in range(d + 1) if ints[j * (d + 1) + k]}
    # held-out verification at primes[2]
    p3 = primes[2]
    f3 = fiber_mod(spec, u_base, ray_class, depth, p3)
    ver = True
    for n in range(depth + 1 - r):
        s = 0
        for j in range(r + 1):
            cj = sum(coeffs.get((j, k), 0) * pow(n, k, p3)
                     for k in range(d + 1)) % p3
            s = (s + cj * f3[n + j]) % p3
        if s != 0:
            ver = False
            break
    log(f"[fiber] ray={ray_class} (r,d)=({r},{d}) "
        f"held-out prime verify={ver}")
    return {'ok': bool(ver), 'ray_class': ray_class, 'u_base': u_base,
            'r': r, 'd': d, 'depth': depth,
            'coeffs': {f'{j},{k}': str(c) for (j, k), c in coeffs.items()},
            'max_coeff_bits': max(abs(c).bit_length()
                                  for c in coeffs.values()),
            'heldout_prime': p3, 'verified': bool(ver)}


def extend_mod(Zs, coeffs, r, d, p, to_n):
    """Extend a fiber sequence mod p to index to_n with a fitted recurrence
    sum_j c_j(n) Z(n+j) = 0.  coeffs keys may be (j,k) or 'j,k' strings.
    Raises on leading-coefficient vanishing (singular index)."""
    cf = {}
    for k, v in coeffs.items():
        jk = tuple(int(t) for t in k.split(',')) if isinstance(k, str) else k
        cf[jk] = int(v)
    Zs = [int(z) % p for z in Zs]
    assert len(Zs) >= r, "need at least r seed values"
    while len(Zs) < to_n + 1:
        n0 = len(Zs) - r
        lead = sum(cf.get((r, k), 0) * pow(n0, k, p)
                   for k in range(d + 1)) % p
        if lead == 0:
            raise ZeroDivisionError(f"leading coeff vanishes at n={n0}")
        s = 0
        for j in range(r):
            cj = sum(cf.get((j, k), 0) * pow(n0, k, p)
                     for k in range(d + 1)) % p
            s = (s + cj * Zs[n0 + j]) % p
        Zs.append((-s) * pow(lead, p - 2, p) % p)
    return Zs


# ---------------------------------------------------- shift-Ore calculus --
# Operators live in the shift algebra Q(n)[S], S: f(n) -> f(n+1), acting on
# fiber sequences as sum_j c_j(n) Z(n+j).  Internal representation: list of
# coefficient polynomials [c_0, .., c_r], each a low-to-high list of
# Fractions; [] is the zero polynomial / zero operator.  Public functions
# take and return the fiber convention (emit_fiber_recurrence-style dicts,
# (coeffs, r, d) tuples, or bare coeffs dicts with (j,k) or 'j,k' keys).

def _pnorm(c):
    c = [Fraction(x) for x in c]
    while c and c[-1] == 0:
        c.pop()
    return c


def _padd(a, b):
    m = max(len(a), len(b))
    return _pnorm([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                   for i in range(m)])


def _pmul(a, b):
    if not a or not b:
        return []
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return _pnorm(out)


def _pshift(a, s):
    """p(n) -> p(n+s)."""
    out, pw = [], [Fraction(1)]
    for ak in a:
        if ak:
            out = _padd(out, [ak * c for c in pw])
        pw = _pmul(pw, [Fraction(s), Fraction(1)])
    return out


def _pdivmod(a, b):
    a, b = _pnorm(a), _pnorm(b)
    assert b, "polynomial division by zero"
    q = [Fraction(0)] * max(0, len(a) - len(b) + 1)
    r = list(a)
    while r and len(r) >= len(b):
        f = r[-1] / b[-1]
        k = len(r) - len(b)
        q[k] = f
        for i in range(len(b)):
            r[k + i] -= f * b[i]
        r = _pnorm(r)
    return _pnorm(q), r


def _pgcd(a, b):
    a, b = _pnorm(a), _pnorm(b)
    while b:
        _, r = _pdivmod(a, b)
        a, b = b, r
    return [c / a[-1] for c in a] if a else []


def _opnorm(A):
    A = [_pnorm(c) for c in A]
    while A and not A[-1]:
        A.pop()
    return A


def _oporder(A):
    return len(A) - 1


def _opsub(A, B):
    m = max(len(A), len(B))
    return _opnorm([_padd(A[i] if i < len(A) else [],
                          [-c for c in (B[i] if i < len(B) else [])])
                    for i in range(m)])


def _opmulpoly(c, A):
    """Left-multiply the operator by a scalar polynomial c(n)."""
    return _opnorm([_pmul(c, aj) for aj in A])


def _opmul(A, B):
    """Operator product A·B (apply B first): S^i·b(n) = b(n+i)·S^i."""
    A, B = _opnorm(A), _opnorm(B)
    if not A or not B:
        return []
    out = [[] for _ in range(len(A) + len(B) - 1)]
    for i, ai in enumerate(A):
        if ai:
            for j, bj in enumerate(B):
                if bj:
                    out[i + j] = _padd(out[i + j], _pmul(ai, _pshift(bj, i)))
    return _opnorm(out)


def _op_primitive(A):
    """Strip the common polynomial factor and rational content; integer
    coefficients, sign-normalized (leading polynomial's leading coeff > 0)."""
    A = _opnorm(A)
    if not A:
        return []
    g = []
    for c in A:
        if c:
            g = _pgcd(g, c)
    if len(g) > 1:
        A = [_pdivmod(c, g)[0] if c else [] for c in A]
    L, G = 1, 0
    for c in A:
        for x in c:
            if x:
                L = L * x.denominator // math.gcd(L, x.denominator)
                G = math.gcd(G, abs(x.numerator))
    s = Fraction(L, max(G, 1))
    A = [[x * s for x in c] for c in A]
    if A[-1][-1] < 0:
        A = [[-x for x in c] for c in A]
    return _opnorm(A)


def _op_prem(A, B):
    """Right pseudo-remainder: some unit c(n) in Q(n) has c·A = Q·B + R with
    order(R) < order(B); returns R (content-stripped — units are free)."""
    A, B = _opnorm(A), _opnorm(B)
    assert B, "pseudo-remainder by the zero operator"
    rb = _oporder(B)
    while A and _oporder(A) >= rb:
        t = _oporder(A) - rb
        killer = [[] for _ in range(t)] + [A[-1]]        # lead_A(n)·S^t
        A = _op_primitive(_opsub(_opmulpoly(_pshift(B[-1], t), A),
                                 _opmul(killer, B)))
    return A


def _op_gcrd(A, B):
    A, B = _op_primitive(A), _op_primitive(B)
    if not A:
        return B
    if not B:
        return A
    if _oporder(A) < _oporder(B):
        A, B = B, A
    while B:
        A, B = B, _op_prem(A, B)
    return A


def _row_scalar_strip(ops):
    """Common rational-content strip across a row of operators (preserves
    any linear identity among them)."""
    L, G = 1, 0
    for op in ops:
        for c in op:
            for x in c:
                if x:
                    L = L * x.denominator // math.gcd(L, x.denominator)
                    G = math.gcd(G, abs(x.numerator))
    s = Fraction(L, max(G, 1))
    return [[[x * s for x in c] for c in op] for op in ops]


def _op_lclm(A, B):
    """LCLM by the extended right Euclidean algorithm: rows keep
    U_i·A + V_i·B = R_i; the zero-remainder row gives LCLM = U·A."""
    A, B = _op_primitive(A), _op_primitive(B)
    assert A and B, "lclm of a zero operator"
    if _oporder(A) < _oporder(B):
        A, B = B, A
    one = [[Fraction(1)]]
    rows = [(A, one, []), (B, [], one)]
    while rows[-1][0]:
        (R0, U0, V0), (R1, U1, V1) = rows[-2], rows[-1]
        W, U, V = R0, U0, V0
        while W and _oporder(W) >= _oporder(R1):
            t = _oporder(W) - _oporder(R1)
            bsh = _pshift(R1[-1], t)
            killer = [[] for _ in range(t)] + [W[-1]]
            W = _opsub(_opmulpoly(bsh, W), _opmul(killer, R1))
            U = _opsub(_opmulpoly(bsh, U), _opmul(killer, U1))
            V = _opsub(_opmulpoly(bsh, V), _opmul(killer, V1))
        W, U, V = _row_scalar_strip([W, U, V])
        rows.append((_opnorm(W), _opnorm(U), _opnorm(V)))
    _, U, _ = rows[-1]
    return _op_primitive(_opmul(U, A))


def _as_op(x):
    """Operator input -> internal list of coefficient polynomials.  Accepts
    an emit_fiber_recurrence-style dict (key 'coeffs'), a (coeffs, r, d)
    tuple, a bare coeffs dict with (j,k) or 'j,k' keys, or the internal
    representation itself."""
    if isinstance(x, list):
        return _opnorm(x)
    if isinstance(x, dict) and 'coeffs' in x:
        coeffs = x['coeffs']
    elif isinstance(x, tuple):
        coeffs = x[0]
    else:
        coeffs = x
    terms = {}
    for k, v in coeffs.items():
        j, kk = (tuple(int(t) for t in k.split(','))
                 if isinstance(k, str) else (int(k[0]), int(k[1])))
        terms.setdefault(j, {})[kk] = Fraction(int(v))
    if not terms:
        return []
    A = []
    for j in range(max(terms) + 1):
        cj = terms.get(j, {})
        A.append([cj.get(k, Fraction(0)) for k in range(max(cj, default=0) + 1)])
    return _opnorm(A)


def _op_dict(A):
    """Internal operator -> fiber-convention dict (primitive integer
    coefficients, emit_fiber_recurrence key style)."""
    A = _op_primitive(A)
    if not A:
        return {'r': -1, 'd': 0, 'coeffs': {}, 'max_coeff_bits': 0}
    coeffs, d = {}, 0
    for j, c in enumerate(A):
        for k, x in enumerate(c):
            if x:
                assert x.denominator == 1
                coeffs[f'{j},{k}'] = str(int(x))
                d = max(d, k)
    return {'r': _oporder(A), 'd': d, 'coeffs': coeffs,
            'max_coeff_bits': max(abs(int(x)).bit_length()
                                  for c in A for x in c if x)}


def op_mul(a, b):
    """Shift-operator product a·b (b applied first), fiber convention."""
    return _op_dict(_opmul(_as_op(a), _as_op(b)))


def right_divides(numer, factor):
    """True iff factor right-divides numer over Q(n) (zero right
    pseudo-remainder — an exact certificate, no floats)."""
    F = _as_op(factor)
    assert F, "right_divides: zero factor"
    return not _op_prem(_as_op(numer), F)


def gcrd(a, b):
    """Primitive greatest common right divisor of two lifted operators —
    the common right factor annihilating any sequence both annihilate."""
    return _op_dict(_op_gcrd(_as_op(a), _as_op(b)))


def lclm(a, b):
    """Primitive least common left multiple (annihilates every sequence
    either operator annihilates); certified by exact right division."""
    A, B = _as_op(a), _as_op(b)
    L = _op_lclm(A, B)
    assert right_divides(L, A) and right_divides(L, B), \
        "lclm certificate failed"
    return _op_dict(L)


def annihilates_mod(op, Zs, p):
    """True iff sum_j c_j(n) Z(n+j) == 0 mod p at every window start n."""
    A = _as_op(op)
    r = _oporder(A)
    if r < 0 or len(Zs) < r + 2:
        return False
    for n0 in range(len(Zs) - r):
        s = 0
        for j in range(r + 1):
            cj, v = A[j], 0
            for coef in reversed(cj):
                v = (v * n0 + coef.numerator
                     * pow(coef.denominator, p - 2, p)) % p
            s = (s + v * (int(Zs[n0 + j]) % p)) % p
        if s:
            return False
    return True


def gcrd_reduce(coeffs, r, d, other=None):
    """GCRD/right-factor reduction of a lifted scalar operator.

    With `other` (a second lifted operator of the same fiber — another
    fit, or a lift from another ray — as an emit_fiber_recurrence-style
    dict, a (coeffs, r, d) tuple, or a bare coeffs dict): returns the
    primitive integer GCRD of the two shift operators, the common right
    factor that annihilates any sequence both inputs annihilate; `ok`
    certifies exact right division of BOTH inputs by the result and
    `reduced` reports a proper order drop.  Without `other`: returns the
    operator content-stripped and sign-normalized (primitive form).

    (r, d) are the operator's order and coefficient degree as emitted;
    the reduction is exact over Q(n) — no fitting, no floats."""
    A = _as_op((coeffs, r, d))
    if other is None:
        out = _op_dict(A)
        out['ok'] = True
        return out
    B = _as_op(other)
    G = _op_gcrd(A, B)
    out = _op_dict(G)
    out['ok'] = bool(G) and right_divides(A, G) and right_divides(B, G)
    out['reduced'] = bool(G) and _oporder(G) < _oporder(_opnorm(A))
    return out

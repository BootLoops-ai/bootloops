"""boxwalk.core — problem spec, dense moment oracle, window primitives.

Objects:
  Z(u) = int_{[0,1]^n} prod_k P_k(x)^{u_k} dx,   P_k in Z[x_1..x_n], u in N^m.
  M(beta;u) = int x^beta prod_k P_k^{u_k} dx  (box moments; Z = M(0)).

Oracle: dense coefficient array Q of prod P_k^{u_k} mod p (sparse-support
convolution), then M(beta) = sum_i Q[i] prod_e 1/(i_e+beta_e+1) via a
separable contraction with the int64 hi/lo overflow split (pilot p0 pattern:
split the accumulating tensor into 16-bit halves so partial sums stay <2^63).

All mod-p arithmetic assumes p < 2^31 (products < 2^62 fit int64).
"""
import json, hashlib, itertools, math
import numpy as np

INT = np.int64


# ---------------------------------------------------------------- spec ----
class Spec:
    """Problem spec. JSON fields:
      name, variables: [str], polynomials: {label: {"e1 e2 .. en": coeff}},
      target_u: {label: int}, order: 'auto' | [labels], options: {}.
    """

    def __init__(self, d):
        self.raw = d
        self.name = d.get('name', 'unnamed')
        self.vars = list(d['variables'])
        self.n = len(self.vars)
        assert 1 <= self.n <= 6, "boxwalk is practical for n <= 6"
        self.polys = {}
        for lab, md in d['polynomials'].items():
            pd = {}
            for ms, c in md.items():
                key = tuple(int(t) for t in ms.split())
                assert len(key) == self.n, f"bad monomial {ms} for n={self.n}"
                if int(c):
                    pd[key] = int(c)
            assert pd, f"empty polynomial {lab}"
            self.polys[lab] = pd
        self.target_u = {k: int(v) for k, v in d.get('target_u', {}).items()
                         if int(v) > 0}
        for k in self.target_u:
            assert k in self.polys, f"target class {k} not in polynomials"
        self.order = d.get('order', 'auto')
        self.options = d.get('options', {})

    @classmethod
    def load(cls, path):
        return cls(json.load(open(path)))

    def degs(self, lab):
        """Per-axis max degree of P_lab."""
        return [max(m[e] for m in self.polys[lab]) for e in range(self.n)]

    def maxdeg(self, lab):
        return max(self.degs(lab))

    def to_dict(self):
        return {'name': self.name, 'variables': self.vars,
                'polynomials': {lab: {' '.join(map(str, m)): c
                                      for m, c in pd.items()}
                                for lab, pd in self.polys.items()},
                'target_u': self.target_u,
                'order': self.order, 'options': self.options}

    def sha(self):
        s = json.dumps(self.to_dict(), sort_keys=True)
        return hashlib.sha256(s.encode()).hexdigest()[:16]


# ------------------------------------------------------------- mod utils --
def modinv(a, p):
    return pow(int(a) % p, p - 2, p)


def primes31(k):
    """k distinct primes descending from 2^31-1 (int64-safe)."""
    from sympy import prevprime
    out, p = [], (1 << 31)
    for _ in range(k):
        p = prevprime(p)
        out.append(int(p))
    return out


def crt_list(residues, moduli):
    """CRT combine; returns (r, M)."""
    r, M = 0, 1
    for a, p in zip(residues, moduli):
        t = (a - r) * pow(M % p, -1, p) % p
        r, M = r + M * t, M * p
    return r % M, M


def ratrec(a, m):
    """Wang rational reconstruction of a mod m -> Fraction or None."""
    from fractions import Fraction
    a %= m
    u0, u1, v0, v1 = m, a, 0, 1
    bound = math.isqrt(m // 2)
    while u1 > bound:
        q = u0 // u1
        u0, u1 = u1, u0 - q * u1
        v0, v1 = v1, v0 - q * v1
    if v1 == 0 or abs(v1) > bound or math.gcd(u1, abs(v1)) != 1:
        return None
    return Fraction(u1, v1) if v1 > 0 else Fraction(-u1, -v1)


# ------------------------------------------------------------ oracle ------
def conv_mod(Q, poly, p):
    """Multiply dense coeff array Q by sparse polynomial dict, mod p."""
    n = Q.ndim
    degs = [max(m[e] for m in poly) for e in range(n)]
    out = np.zeros(tuple(s + d for s, d in zip(Q.shape, degs)), dtype=INT)
    for mo, c in poly.items():
        cc = int(c) % p
        if cc == 0:
            continue
        sl = tuple(slice(mo[e], mo[e] + Q.shape[e]) for e in range(n))
        out[sl] = (out[sl] + cc * Q) % p
    return out


def dense_Q(spec, u, p):
    """Dense coefficient array of prod P_k^{u_k} mod p."""
    Q = np.ones((1,) * spec.n, dtype=INT)
    for lab, e in u.items():
        for _ in range(int(e)):
            Q = conv_mod(Q, spec.polys[lab], p)
    return Q


def _mod_tensordot0(T, W, p):
    """tensordot(T, W, axes=([0],[0])) mod p, int64-safe hi/lo split."""
    T = T % p
    Thi, Tlo = T >> 16, T & 0xFFFF
    hi = np.tensordot(Thi, W, axes=([0], [0])) % p
    lo = np.tensordot(Tlo, W, axes=([0], [0])) % p
    return ((hi << 16) + lo) % p


def window_moments(Q, bmax, p):
    """All M(beta) for beta in [0,bmax]^n from dense Q mod p."""
    n = Q.ndim
    mmax = max(Q.shape)
    inv = np.array([0] + [pow(j, p - 2, p) for j in range(1, mmax + bmax + 1)],
                   dtype=INT)
    T = Q % p
    for _ in range(n):
        m = T.shape[0]
        idx = np.arange(m)[:, None] + np.arange(bmax + 1)[None, :] + 1
        W = inv[idx]
        T = _mod_tensordot0(T, W, p)     # contracted axis result goes last
    return T


def moment_window(spec, u, bmax, p):
    return window_moments(dense_Q(spec, u, p), bmax, p)


def Z_mod(spec, u, p):
    """Exact Z(u) mod p via the dense oracle."""
    return int(moment_window(spec, u, 0, p).ravel()[0])


# ------------------------------------------------------ window primitives -
def init_window_batch(side, n, pvec):
    """M(beta;0) = prod 1/(beta_e+1), batched over primes: (K,)+(side,)*n."""
    K = len(pvec)
    M = np.ones((K,) + (side,) * n, dtype=INT)
    for k, p in enumerate(pvec):
        p = int(p)
        inv = np.array([pow(j + 1, p - 2, p) for j in range(side)], dtype=INT)
        Mk = np.ones((side,) * n, dtype=INT)
        for ax in range(n):
            sh = [1] * n
            sh[ax] = side
            Mk = (Mk * inv.reshape(sh)) % p
        M[k] = Mk
    return M


def raise_batch(M, poly, pvec, valid):
    """One contiguity raise u -> u+e_k, prime-batched.
    M: (K,)+(side,)*n with entries valid on [0,valid_e); returns (Mnew,
    newvalid) with newvalid_e = valid_e - deg_e(P_k); rim set to -1.
    M'(beta) = sum_gamma P_k[gamma] * M(beta+gamma)."""
    n = M.ndim - 1
    K = M.shape[0]
    degs = [max(m[e] for m in poly) for e in range(n)]
    newvalid = [valid[e] - degs[e] for e in range(n)]
    if min(newvalid) < 1:
        raise RuntimeError(f"window exhausted: valid={valid} degs={degs}")
    pv = pvec.reshape((K,) + (1,) * n)
    acc = np.zeros((K,) + tuple(newvalid), dtype=INT)
    for mo, c in poly.items():
        cv = np.array([int(c) % int(p) for p in pvec],
                      dtype=INT).reshape((K,) + (1,) * n)
        sl = (slice(None),) + tuple(slice(mo[e], mo[e] + newvalid[e])
                                    for e in range(n))
        acc = (acc + cv * M[sl]) % pv
    out = np.full_like(M, -1)
    out[(slice(None),) + tuple(slice(0, newvalid[e]) for e in range(n))] = acc
    return out, newvalid

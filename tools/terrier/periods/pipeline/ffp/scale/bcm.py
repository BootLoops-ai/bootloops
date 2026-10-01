# BCM finite-hypergeometric a_p(z) engine (Gauss-sum definition over F_q, complex-exact:
# H_q is a rational integer; float residual certified < 1e-4 before rounding).
# Requires (q-1)*a_i in Z (choose p = 1 mod lcm(denoms) for F_p pass; q = p^2 pass works
# for denoms | 24 at any odd p coprime to them). beta = (1,...,1) (MUM ops).
import numpy as np
from fractions import Fraction

def _factor(n):
    fs, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            fs[d] = fs.get(d, 0) + 1; n //= d
        d += 1
    if n > 1:
        fs[n] = fs.get(n, 0) + 1
    return fs

class GaussField:
    """F_q Gauss-sum table, q = p^deg (deg 1 or 2). G[r] = sum_x w^r(x) psi(x),
    w(g^k) = zeta_{q-1}^k, psi = e(Tr(x)/p). Also dlog of F_p^* elements."""
    def __init__(self, p, deg):
        self.p, self.deg, self.q = p, deg, p ** deg
        N = self.q - 1
        if deg == 1:
            g = self._gen_fp(p)
            tr = np.empty(N, dtype=np.int64); pos_fp = np.zeros(p, dtype=np.int64)
            x = 1
            for k in range(N):
                tr[k] = x; pos_fp[x] = k; x = x * g % p
        else:
            n = 2
            while pow(n, (p - 1) // 2, p) != p - 1:
                n += 1
            g = self._gen_fp2(p, n)
            tr = np.empty(N, dtype=np.int64); pos_fp = np.zeros(p, dtype=np.int64)
            a, b = 1, 0
            for k in range(N):
                tr[k] = 2 * a % p
                if b == 0:
                    pos_fp[a] = k
                a, b = (a * g[0] + n * b * g[1]) % p, (a * g[1] + b * g[0]) % p
        self.pos_fp = pos_fp
        v = np.exp(2j * np.pi * tr / p)
        self.G = np.fft.ifft(v) * N          # G[r] = sum_k zeta^{rk} v[k]
        assert abs(self.G[0] - (-1)) < 1e-8, "g(triv) != -1"
        gen_mag = np.abs(self.G[1:])
        assert np.max(np.abs(gen_mag - np.sqrt(self.q))) < 1e-6, "Gauss magnitude"

    @staticmethod
    def _gen_fp(p):
        fs = list(_factor(p - 1))
        for g in range(2, p):
            if all(pow(g, (p - 1) // f, p) != 1 for f in fs):
                return g
        raise ValueError

    @staticmethod
    def _gen_fp2(p, n):
        N = p * p - 1
        fs = list(_factor(N))
        def powe(el, e):
            ra, rb = 1, 0; a, b = el
            while e:
                if e & 1:
                    ra, rb = (ra * a + n * rb * b) % p, (ra * b + rb * a) % p
                a, b = (a * a + n * b * b) % p, (2 * a * b) % p
                e >>= 1
            return ra, rb
        for aa in range(1, p):
            for bb in range(1, p):
                if all(powe((aa, bb), N // f) != (1, 0) for f in fs):
                    return (aa, bb)
        raise ValueError

def h_coeffs(gf, alpha):
    """C[m] for beta=(1,)*len(alpha); needs (q-1)*a_i integral."""
    N = gf.q - 1
    A = []
    for a in alpha:
        fr = Fraction(a)
        assert N % fr.denominator == 0, (gf.q, fr, 'q-1 not divisible by denom')
        A.append(N * fr.numerator // fr.denominator % N)
    m = np.arange(N)
    C = np.ones(N, dtype=np.complex128)
    for Ai in A:
        C *= gf.G[(m + Ai) % N] / gf.G[Ai]
    C *= (gf.G[(N - m) % N] / gf.G[0]) ** len(alpha)
    return C

def h_all_fp(gf, C, conj=False):
    """H(t) for all t in F_p^* (index by t). Returns complex array len p (0 unused)."""
    N, p = gf.q - 1, gf.p
    zeta = np.exp((-2j if conj else 2j) * np.pi * np.arange(N) / N)
    m = np.arange(N)
    out = np.zeros(p, dtype=np.complex128)
    for t in range(1, p):
        kt = int(gf.pos_fp[t])
        out[t] = np.sum(C * zeta[(m * kt) % N])
    return out / (1 - gf.q)

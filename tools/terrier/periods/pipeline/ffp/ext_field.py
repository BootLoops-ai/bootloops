# F_{p^m} arithmetic on numpy arrays of coefficient vectors (batched), for mu_lambda measures.
import numpy as np
import sympy

def find_irred(p, m):
    """Monic irreducible of degree m over F_p; returns low-part coeffs of x^m (i.e. x^m = low)."""
    x = sympy.symbols('x')
    for c in range(1, 3 * p):
        a, b = c % p, (c // p) % p
        f = sympy.Poly(x**m - b * x - a, x, modulus=p)
        if f.is_irreducible:
            return np.array([a % p, b % p] + [0] * (m - 2), dtype=np.int64)
    raise ValueError("no irreducible found")

class Fpm:
    """Elements: int64 arrays shaped (..., m), coeffs of 1..x^{m-1} mod f, x^m = low (deg<=1)."""
    def __init__(self, p, m):
        self.p, self.m = p, m
        low = find_irred(p, m)
        red = np.zeros((m - 1, m), dtype=np.int64)
        cur = low.copy()
        red[0] = cur
        for j in range(1, m - 1):
            nxt = np.roll(cur, 1); hi = int(nxt[0]); nxt[0] = 0
            nxt = (nxt + hi * low) % p
            red[j] = nxt; cur = nxt
        self.low, self.red = low, red

    def mul(self, u, v):
        p, m = self.p, self.m
        prod = np.zeros(np.broadcast_shapes(u.shape[:-1], v.shape[:-1]) + (2 * m - 1,),
                        dtype=np.int64)
        for i in range(m):
            for j in range(m):
                prod[..., i + j] += u[..., i] * v[..., j]   # <= 25 p^2 * ... safe in int64
        prod %= p
        out = prod[..., :m].copy()
        for j in range(m - 1):
            out = out + prod[..., m + j:m + j + 1] * self.red[j]
        return out % p

    def inv(self, u):
        e = self.p ** self.m - 2
        r = np.zeros_like(u); r[..., 0] = 1
        b = u % self.p
        while e:
            if e & 1:
                r = self.mul(r, b)
            b = self.mul(b, b)
            e >>= 1
        return r

    def trace_vec(self):
        p, m = self.p, self.m
        tv = np.zeros(m, dtype=np.int64)
        # Tr(x^j): power sums of roots of f. Newton's identities on f = x^m - b x - a.
        # f coeffs: e_k = 0 except e_{m-1} = -(-b)?? Use direct: companion-matrix power traces.
        C = np.zeros((m, m), dtype=np.int64)
        for i in range(1, m):
            C[i, i - 1] = 1
        C[:, m - 1] = self.low  # x * x^{m-1} = x^m = low
        M = np.eye(m, dtype=np.int64)
        for j in range(m):
            tv[j] = np.trace(M) % p
            M = (M @ C) % p
        return tv

    def trace(self, u):
        if not hasattr(self, '_trv'):
            self._trv = self.trace_vec()
        return (u * self._trv).sum(axis=-1) % self.p

def batch_inv(K, y, blk=512):
    """Montgomery blocked inversion of nonzero elements, shape (n, m)."""
    n = y.shape[0]
    pad = (-n) % blk
    if pad:
        onep = np.zeros((pad, K.m), dtype=np.int64); onep[:, 0] = 1
        y = np.vstack([y, onep])
    rows = y.reshape(-1, blk, K.m)
    pref = np.empty_like(rows)
    pref[:, 0] = rows[:, 0]
    for j in range(1, blk):
        pref[:, j] = K.mul(pref[:, j - 1], rows[:, j])
    tot_inv = K.inv(pref[:, blk - 1])          # Fermat on n/blk elements only
    out = np.empty_like(rows)
    acc = tot_inv
    for j in range(blk - 1, 0, -1):
        out[:, j] = K.mul(acc, pref[:, j - 1])
        acc = K.mul(acc, rows[:, j])
    out[:, 0] = acc
    res = out.reshape(-1, K.m)
    return res[:y.shape[0] - pad] if pad else res

def mu_row(p, lam, tr_val=1, chunk=1 << 20):
    """F(t) = #{y in F_{p^lam}^*: Tr y = tr_val, Tr 1/y = t}: enumerate the affine slice
    (p^{lam-1} elements) parametrized by (c_1..c_{lam-1}); c_0 solved from the trace."""
    K = Fpm(p, lam)
    tv = K.trace_vec()
    assert tv[0] % p != 0, "Tr(1) = 0 mod p — bad degree/prime combo"
    inv_t0 = pow(int(tv[0]), p - 2, p)
    nslice = p ** (lam - 1)
    cnt = np.zeros(p, dtype=np.int64)
    for lo in range(0, nslice, chunk):
        idx = np.arange(lo, min(lo + chunk, nslice), dtype=np.int64)
        y = np.zeros((len(idx), lam), dtype=np.int64)
        t = idx.copy()
        for j in range(1, lam):
            y[:, j] = t % p; t //= p
        s = (y[:, 1:] * tv[1:]).sum(axis=1) % p
        y[:, 0] = ((tr_val - s) * inv_t0) % p
        nz = np.any(y != 0, axis=1)
        y = y[nz]
        tri = K.trace(batch_inv(K, y))
        cnt += np.bincount(tri, minlength=p)
    return cnt

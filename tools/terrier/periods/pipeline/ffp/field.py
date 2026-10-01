# F_q index arithmetic, q = p or p^2. Elements are indices 0..q-1.
# q=p: index = value. q=p^2: index a*p+b <-> a + b*w, w^2 = n (n = least nonresidue).
import numpy as np

def inv_table_fp(p):
    inv = np.zeros(p, dtype=np.int64)
    inv[1] = 1
    for i in range(2, p):
        inv[i] = (-(p // i) * inv[p % i]) % p
    return inv

def least_nonresidue(p):
    sq = set((i * i) % p for i in range(1, p))
    for n in range(2, p):
        if n not in sq:
            return n
    raise ValueError

class Fq:
    def __init__(self, p, deg):
        assert deg in (1, 2)
        self.p, self.deg, self.q = p, deg, p ** deg
        self.invp = inv_table_fp(p)
        q = self.q
        idx = np.arange(q, dtype=np.int64)
        if deg == 1:
            self.A, self.B = idx, None
            self.inv = np.zeros(q, dtype=np.int64)
            self.inv[1:] = self.invp[1:]
        else:
            self.n = least_nonresidue(p)
            a, b = idx // p, idx % p
            self.A, self.B = a, b
            # inverse: (a+bw)^-1 = (a - bw)/(a^2 - n b^2)
            norm = (a * a - self.n * b * b) % p
            invn = np.zeros(q, dtype=np.int64)
            invn[1:] = self.invp[norm[1:]]
            self.inv = ((a * invn) % p) * p + ((-b * invn) % p)
            self.inv[0] = 0

    def embed_fp(self, v):
        """index of the prime-field element v (int)."""
        v %= self.p
        return v if self.deg == 1 else v * self.p

    def mul_scalar_all(self, s):
        """s * t for all t in F_q; s an index. Returns index array."""
        p = self.p
        if self.deg == 1:
            return (s * self.A) % p
        s1, s2 = s // p, s % p
        c1 = (s1 * self.A + self.n * s2 * self.B) % p
        c2 = (s1 * self.B + s2 * self.A) % p
        return c1 * p + c2

    def add_scalar(self, arr, c):
        """arr (+) c elementwise, c an index."""
        p = self.p
        if self.deg == 1:
            return (arr + c) % p
        c1, c2 = c // p, c % p
        return ((arr // p + c1) % p) * p + ((arr % p + c2) % p)

    def sub(self, i, j):
        p = self.p
        if self.deg == 1:
            return (i - j) % p
        return ((i // p - j // p) % p) * p + ((i % p - j % p) % p)

    def one(self):
        return 1 if self.deg == 1 else self.p  # 1 + 0w -> a=1,b=0 -> idx = p

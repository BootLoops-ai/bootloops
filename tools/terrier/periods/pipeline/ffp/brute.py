# Brute-force validators for small q.
import numpy as np
from field import Fq
from hv_count import torus_counts, twisted12_counts

def brute_torus(p, deg):
    K = Fq(p, deg)
    q, one = K.q, K.one()
    def add(i, j):
        if deg == 1:
            return (i + j) % p
        return ((i // p + j // p) % p) * p + ((i % p + j % p) % p)
    def neg(i):
        return K.sub(0, i)
    N = np.zeros(q, dtype=np.int64)
    units = [x for x in range(1, q)]
    for x1 in units:
        for x2 in units:
            s12 = add(x1, x2)
            i12 = add(K.inv[x1], K.inv[x2])
            for x3 in units:
                s3 = add(s12, x3)
                i3 = add(i12, K.inv[x3])
                for x4 in units:
                    x5 = K.sub(one, add(s3, x4))
                    if x5 == 0:
                        continue
                    psi = add(add(i3, K.inv[x4]), K.inv[x5])
                    N[psi] += 1
    return N

def brute_twisted12(p):
    Kp, K2 = Fq(p, 1), Fq(p, 2)
    N = np.zeros(p, dtype=np.int64)
    for y in range(1, p * p):
        tr = (2 * (y // p)) % p
        tri = (2 * (K2.inv[y] // p)) % p
        for x3 in range(1, p):
            for x4 in range(1, p):
                x5 = (1 - tr - x3 - x4) % p
                if x5 == 0:
                    continue
                psi = (tri + Kp.inv[x3] + Kp.inv[x4] + Kp.inv[x5]) % p
                N[psi] += 1
    return N

if __name__ == '__main__':
    for (p, deg) in [(7, 1), (11, 1), (13, 1), (5, 2)]:
        K, F = torus_counts(p, deg)
        B = brute_torus(p, deg)
        ok = np.array_equal(F, B)
        print(f"torus p={p} deg={deg}: {'OK' if ok else 'MISMATCH'}",
              F[:6].tolist(), B[:6].tolist() if not ok else '')
        assert ok
    for p in [5, 7, 11]:
        F = twisted12_counts(p)
        B = brute_twisted12(p)
        ok = np.array_equal(F, B)
        print(f"twist12 p={p}: {'OK' if ok else 'MISMATCH'}", '' if ok else (F.tolist(), B.tolist()))
        assert ok
    print("ALL BRUTE CHECKS PASS")

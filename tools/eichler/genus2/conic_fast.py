# eichler.genus2 — conic pipeline in exact (p,q) pair arithmetic over K=Q(sqrt409). No sympy simplify.
import json, itertools
from fractions import Fraction as Fr

D409 = 409
class K:
    __slots__ = ('p', 'q')
    def __init__(s, p=0, q=0): s.p = Fr(p); s.q = Fr(q)
    def __add__(a, b): return K(a.p+b.p, a.q+b.q)
    def __sub__(a, b): return K(a.p-b.p, a.q-b.q)
    def __neg__(a): return K(-a.p, -a.q)
    def __mul__(a, b):
        if isinstance(b, (int, Fr)): return K(a.p*b, a.q*b)
        return K(a.p*b.p + D409*a.q*b.q, a.p*b.q + a.q*b.p)
    __rmul__ = __mul__
    def inv(a):
        n = a.p*a.p - D409*a.q*a.q
        return K(a.p/n, -a.q/n)
    def __truediv__(a, b):
        if isinstance(b, (int, Fr)): return K(a.p/b, a.q/b)
        return a*b.inv()
    def iszero(a): return a.p == 0 and a.q == 0
    def __repr__(a): return f'({a.p} + {a.q}*s409)'
    def conj(a): return K(a.p, -a.q)

# exact invariants
j1 = K(Fr(556889505536, 225), Fr(27279625472, 225))
j2 = K(Fr(31137174272, 225), Fr(1532443904, 225))
j3 = K(Fr(2646125216, 75), Fr(130109792, 75))
I2 = K(1)
I4 = j2/j1
I6 = j3/j1
I10 = K(1)/j1
# Mestre x,y,z (I2=1 simplifies)
x = (K(1) + I4*20)* Fr(8, 225)
y = (K(1) + I4*80 - I6*600) * Fr(16, 3375)
z = (K(Fr(-10800000))*I10 + K(-9) - I4*700 + I6*3600 + I4*I4*12400 - I4*I6*48000) * Fr(-64, 253125)
print('x =', x); print('y =', y); print('z =', z)
L = [[x + y*6,        x*x*6 + y*2,   z*2],
     [x*x*6 + y*2,    z*2,           x*x*x*9 + x*y*4 + y*y*6],
     [z*2,            x*x*x*9 + x*y*4 + y*y*6, x*x*y*6 + y*y*2 + x*z*3]]
# congruence diagonalization
P = [[K(1) if i==j else K(0) for j in range(3)] for i in range(3)]
A = [row[:] for row in L]
def matmul(M, N):
    return [[sum((M[i][k]*N[k][j] for k in range(3)), K(0)) for j in range(3)] for i in range(3)]
def transpose(M): return [[M[j][i] for j in range(3)] for i in range(3)]
for k in range(3):
    if A[k][k].iszero():
        done = False
        for l in range(k+1, 3):
            if not A[l][l].iszero():
                T = [[K(1) if i==j else K(0) for j in range(3)] for i in range(3)]
                T[k][k] = K(0); T[l][l] = K(0); T[k][l] = K(1); T[l][k] = K(1)
                A = matmul(matmul(T, A), transpose(T)); P = matmul(T, P); done = True; break
        if not done:
            for l in range(k+1, 3):
                if not A[k][l].iszero():
                    T = [[K(1) if i==j else K(0) for j in range(3)] for i in range(3)]
                    T[k][l] = K(1)
                    A = matmul(matmul(T, A), transpose(T)); P = matmul(T, P); break
    piv = A[k][k]
    T = [[K(1) if i==j else K(0) for j in range(3)] for i in range(3)]
    for l in range(k+1, 3):
        T[l][k] = -(A[l][k]/piv)
    A = matmul(matmul(T, A), transpose(T)); P = matmul(T, P)
d = [A[i][i] for i in range(3)]
print('diagonal:')
for i in range(3): print(' d%d =' % i, d[i])
# real-embedding sign check (both embeddings): sigma(sqrt409) = ±20.2237...
import math
s = math.sqrt(409)
for emb, tag in ((s, '+'), (-s, '-')):
    signs = [1 if float(t.p) + float(t.q)*emb > 0 else -1 for t in d]
    print(f'embedding {tag}: signs', signs, '(indefinite needed for real points)')

def is_square_in_K(e):
    if e.iszero(): return K(0)
    p, q = e.p, e.q
    if q == 0:
        if p < 0: return None
        r = p
        # sqrt of Fraction rational?
        import math as m
        def frsqrt(fr):
            if fr < 0: return None
            n, dd = fr.numerator, fr.denominator
            rn = m.isqrt(n); rd = m.isqrt(dd)
            if rn*rn == n and rd*rd == dd: return Fr(rn, rd)
            return None
        rr = frsqrt(p)
        if rr is not None: return K(rr, 0)
        rr = frsqrt(p/D409)
        if rr is not None: return K(0, rr)
        return None
    # (a+b s)^2 = e: a^2+409b^2 = p, 2ab = q: a^2 solves t^2 - p t + 409 q^2/4 = 0
    disc = p*p - D409*q*q
    import math as m
    def frsqrt(fr):
        if fr < 0: return None
        n, dd = fr.numerator, fr.denominator
        rn = m.isqrt(n); rd = m.isqrt(dd)
        if rn*rn == n and rd*rd == dd: return Fr(rn, rd)
        return None
    rd = frsqrt(disc)
    if rd is None: return None
    for t0 in ((p+rd)/2, (p-rd)/2):
        a = frsqrt(t0)
        if a is not None and a != 0:
            b = q/(2*a)
            return K(a, b)
    return None

found = None
R = 60
for v0 in range(-R, R+1):
    for w0 in range(-R, R+1):
        if v0 == 0 and w0 == 0: continue
        rhs = -(d[1]*(v0*v0) + d[2]*(w0*w0))/d[0]
        r = is_square_in_K(rhs)
        if r is not None:
            found = (r, v0, w0)
            break
    if found: break
if found:
    u0, v0, w0 = found
    print('DIAG POINT:', u0, v0, w0)
    # back to original coords: A = P L P^T => point uL = P^T [u0,v0,w0]
    PT = transpose(P)
    vec = [sum((PT[i][k]*[u0, K(v0), K(w0)][k] for k in range(3)), K(0)) for i in range(3)]
    chk = sum((vec[i]*L[i][j]*vec[j] for i in range(3) for j in range(3)), K(0))
    print('back-check (0?):', chk)
    json.dump({'point': [[str(t.p), str(t.q)] for t in vec],
               'diag': [[str(t.p), str(t.q)] for t in d]}, open('conic_point.json', 'w'), indent=1)
else:
    print('NO K-POINT in box', R, '-> Mestre obstruction candidate; recording diagonal + delta')
    delta = -(d[2]/d[0])
    json.dump({'diag': [[str(t.p), str(t.q)] for t in d],
               'delta_ext': [str(delta.p), str(delta.q)]}, open('conic_obstruction.json', 'w'), indent=1)

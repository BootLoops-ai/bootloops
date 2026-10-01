"""Stage-3 coefficient arithmetic.
  GT    : F0 = Q(i)(t), elements (A + i B)/D with A, B, D fmpq_poly in t (D monic, content-normalized).
  KF    : simple extension F0[y]/(f(y)) over F0 = Q(i)(t) + all Laurent coefficients of
          N(y)/(lead prod f_i^m_i) at a root of f_k.
  Points: fibration points = roots in x of alphabet letters: ('Q', GT value) for letters linear in x, ('A', letter, slot) for degree >= 2.
  TRing : the tensor ring F0[y_g : g generators]/(slot relations) in which all stage-3 coefficients live; generators g = (letter, slot);
          slot relations: P_0 = f/lc (slot 0), P_s(y_s) = (P_{s-1}(y_s) - P_{s-1}(y_{s-1}))/(y_s - y_{s-1}) (distinct roots of the same letter).
          Inverses are only ever needed for (i) single-letter elements (extended Euclid in KF) and (ii) point differences a - p and points
          themselves, done by the closed forms documented in inv_point / inv_diff."""
import sys
import slice_prov
from fractions import Fraction as Fr
from flint import fmpq, fmpq_poly, fmpq_mpoly_ctx
from math import comb

_P1 = fmpq_poly([1]); _P0 = fmpq_poly([])
def _fp(x):
    if isinstance(x, fmpq_poly): return x
    if isinstance(x, Fr): return fmpq_poly([fmpq(x.numerator, x.denominator)])
    if isinstance(x, int): return fmpq_poly([x])
    if isinstance(x, fmpq): return fmpq_poly([x])
    raise TypeError(type(x))

class GT:
    """(A + i B)/D in Q(i)(t)"""
    __slots__ = ("A", "B", "D")
    def __init__(self, A=0, B=0, D=1, _norm=True):
        if isinstance(A, GT): self.A, self.B, self.D = A.A, A.B, A.D; return
        A = _fp(A); B = _fp(B); D = _fp(D)
        if _norm:
            if D.degree() < 0: raise ZeroDivisionError("GT zero denominator")
            if A.degree() < 0 and B.degree() < 0: D = _P1
            else:
                g = A.gcd(B).gcd(D)
                if g.degree() > 0: A = divmod(A, g)[0]; B = divmod(B, g)[0]; D = divmod(D, g)[0]
                lc = D.coeffs()[-1]
                if lc != 1: A = A / lc; B = B / lc; D = D / lc
        self.A, self.B, self.D = A, B, D
    @staticmethod
    def t(): return GT(fmpq_poly([0, 1]))
    @staticmethod
    def i(): return GT(0, 1)
    def is_zero(self): return self.A.degree() < 0 and self.B.degree() < 0
    def is_real(self): return self.B.degree() < 0
    def __bool__(self): return not self.is_zero()
    def __eq__(self, o):
        if not isinstance(o, GT): o = GT(o)
        return self.A == o.A and self.B == o.B and self.D == o.D
    def __hash__(self): return hash((str(self.A), str(self.B), str(self.D)))
    def __neg__(self): return GT(-self.A, -self.B, self.D, _norm=False)
    def conj(self): return GT(self.A, -self.B, self.D, _norm=False)
    def __add__(self, o):
        if not isinstance(o, GT): o = GT(o)
        if o.is_zero(): return self
        if self.is_zero(): return o
        if self.D == o.D: return GT(self.A + o.A, self.B + o.B, self.D)
        return GT(self.A * o.D + o.A * self.D, self.B * o.D + o.B * self.D, self.D * o.D)
    __radd__ = __add__
    def __sub__(self, o):
        if not isinstance(o, GT): o = GT(o)
        return self + (-o)
    def __rsub__(self, o): return GT(o) - self
    def __mul__(self, o):
        if not isinstance(o, GT):
            if isinstance(o, (int, Fr)):
                if o == 0: return GT()
                q = _fp(Fr(o)); return GT(self.A * q, self.B * q, self.D, _norm=False) if Fr(o).denominator == 1 else GT(self.A * q, self.B * q, self.D)
            o = GT(o)
        if self.is_zero() or o.is_zero(): return GT()
        if o.B.degree() < 0: return GT(self.A * o.A, self.B * o.A, self.D * o.D)
        if self.B.degree() < 0: return GT(self.A * o.A, self.A * o.B, self.D * o.D)
        return GT(self.A * o.A - self.B * o.B, self.A * o.B + self.B * o.A, self.D * o.D)
    __rmul__ = __mul__
    def inv(self):
        if self.is_zero(): raise ZeroDivisionError("GT inverse of 0")
        n = self.A * self.A + self.B * self.B
        return GT(self.A * self.D, -self.B * self.D, n)
    def __truediv__(self, o):
        if not isinstance(o, GT): o = GT(o)
        return self * o.inv()
    def __rtruediv__(self, o): return GT(o) * self.inv()
    def __pow__(self, e):
        if e < 0: return self.inv() ** (-e)
        r = GT(1); b = self
        while e:
            if e & 1: r = r * b
            b = b * b; e >>= 1
        return r
    def degs(self): return (max(self.A.degree(), self.B.degree(), 0), max(self.D.degree(), 0))
    def __call__(self, tv):
        """exact value at rational t -> (Fraction re, Fraction im)"""
        q = fmpq(Fr(tv).numerator, Fr(tv).denominator); d = self.D(q)
        if d == 0: raise ZeroDivisionError("pole at t=%s" % tv)
        a = self.A(q) / d; b = self.B(q) / d
        return Fr(int(a.p), int(a.q)), Fr(int(b.p), int(b.q))
    def mp(self, tv, mpmod):
        re, im = self(tv); return mpmod.mpc(mpmod.mpf(re.numerator) / re.denominator, mpmod.mpf(im.numerator) / im.denominator)
    def __repr__(self): return "GT((%s + i*(%s))/(%s))" % (str(self.A).replace("x", "t"), str(self.B).replace("x", "t"), str(self.D).replace("x", "t"))
    def ser(self): return [[[int(c.p), int(c.q)] for c in P.coeffs()] for P in (self.A, self.B, self.D)]
    @staticmethod
    def deser(s):
        A, B, D = [fmpq_poly([fmpq(p, q) for p, q in lst]) for lst in s]; return GT(A, B, D)
    @staticmethod
    def from_tj_poly(P, jt, jj):
        """x-free fmpq_mpoly in (..., t, j) reduced mod j^2+1 -> GT"""
        d = P.to_dict(); A = {}; B = {}
        for e, c in d.items():
            assert all(int(z) == 0 for k, z in enumerate(e) if k not in (jt, jj)), "polynomial is not x-free"
            k = int(e[jj]); assert k <= 1
            (A if k == 0 else B)[int(e[jt])] = c
        def up(dd):
            if not dd: return _P0
            m = max(dd); return fmpq_poly([dd.get(i, fmpq(0)) for i in range(m + 1)])
        return GT(up(A), up(B))

def xcoeffs_GT(P, jx, jt, jj):
    """fmpq_mpoly in (x's, t, j) -> list of GT: coefficients of x^k (other x's must be absent)"""
    d = P.to_dict(); byk = {}
    for e, c in d.items():
        assert all(int(z) == 0 for k, z in enumerate(e) if k not in (jx, jt, jj)), "unexpected variable in polynomial"
        kk = int(e[jj]); assert kk <= 1
        slot = byk.setdefault(int(e[jx]), ({}, {})); slot[kk][int(e[jt])] = c
    out = []
    for k in range(max(byk) + 1 if byk else 0):
        A, B = byk.get(k, ({}, {}))
        def up(dd):
            if not dd: return _P0
            m = max(dd); return fmpq_poly([dd.get(i, fmpq(0)) for i in range(m + 1)])
        out.append(GT(up(A), up(B)))
    return out

# ---------------------------------------------------------------- polynomials over F0 (lists of GT, increasing degree)
def ptrim(p):
    p = list(p)
    while p and p[-1].is_zero(): p.pop()
    return p
def pmul(a, b):
    if not a or not b: return []
    out = [GT() for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        if x.is_zero(): continue
        for j, y in enumerate(b):
            if y.is_zero(): continue
            out[i + j] = out[i + j] + x * y
    return ptrim(out)
def padd(a, b):
    n = max(len(a), len(b)); return ptrim([(a[i] if i < len(a) else GT()) + (b[i] if i < len(b) else GT()) for i in range(n)])
def psub(a, b): return padd(a, [-x for x in b])
def pdivmod(a, b):
    a = ptrim(a); b = ptrim(b); assert b, "division by zero polynomial"
    q = [GT() for _ in range(max(len(a) - len(b) + 1, 1))]; r = list(a); ilc = b[-1].inv()
    while len(ptrim(r)) >= len(b):
        r = ptrim(r); k = len(r) - len(b); c = r[-1] * ilc; q[k] = q[k] + c
        for i, bc in enumerate(b): r[i + k] = r[i + k] - c * bc
        r = ptrim(r)
    return ptrim(q), ptrim(r)
def ppow(a, e):
    r = [GT(1)]
    for _ in range(e): r = pmul(r, a)
    return r
def peval(p, x):
    acc = GT()
    for c in reversed(p): acc = acc * x + c
    return acc

class KF:
    """F0[y]/(f), f = list of GT (increasing degree), deg >= 1"""
    def __init__(self, fcoeffs):
        self.f = ptrim([GT(c) for c in fcoeffs]); self.d = len(self.f) - 1; assert self.d >= 1
        self.lc = self.f[-1]; self._T = [tuple(GT(1) if i == 0 else GT() for i in range(self.d))]
    def zero(self): return tuple(GT() for _ in range(self.d))
    def one(self): return self._T[0]
    def const(self, q): v = [GT()] * self.d; v[0] = GT(q); return tuple(v)
    def gen(self): return self.tpow(1) if self.d > 1 else self.tpow(1)
    def tpow(self, n):
        while len(self._T) <= n:
            prev = [GT()] + list(self._T[-1])
            if len(prev) > self.d:
                top = prev[self.d]; prev = prev[:self.d]
                if not top.is_zero():
                    for i in range(self.d): prev[i] = prev[i] - top * self.f[i] / self.lc
            self._T.append(tuple(prev))
        return self._T[n]
    def add(self, a, b): return tuple(x + y for x, y in zip(a, b))
    def sub(self, a, b): return tuple(x - y for x, y in zip(a, b))
    def scal(self, q, a): q = GT(q); return tuple(q * x for x in a)
    def is_zero(self, a): return all(x.is_zero() for x in a)
    def mul(self, a, b):
        prod = [GT()] * (2 * self.d - 1)
        for i, x in enumerate(a):
            if x.is_zero(): continue
            for j, y in enumerate(b):
                if y.is_zero(): continue
                prod[i + j] = prod[i + j] + x * y
        out = [GT()] * self.d
        for n, c in enumerate(prod):
            if c.is_zero(): continue
            if n < self.d: out[n] = out[n] + c
            else:
                tn = self.tpow(n)
                for i in range(self.d): out[i] = out[i] + c * tn[i]
        return tuple(out)
    def from_ypoly(self, coeffs):
        out = [GT()] * self.d
        for n, c in enumerate(coeffs):
            c = GT(c)
            if c.is_zero(): continue
            tn = self.tpow(n)
            for i in range(self.d): out[i] = out[i] + c * tn[i]
        return tuple(out)
    def inv(self, a):
        if self.d == 1: return (a[0].inv(),)
        r0, r1 = ptrim(self.f), ptrim(list(a)); s0, s1 = [], [GT(1)]
        while r1:
            q, r = pdivmod(r0, r1); r0, r1 = r1, r; s0, s1 = s1, psub(s0, pmul(q, s1))
        assert len(r0) == 1, "element not invertible in KF (f reducible over Q(i)(t)? or zero divisor)"
        c = r0[0].inv(); return self.from_ypoly([x * c for x in s0])
    def pow(self, a, e):
        if e < 0: return self.pow(self.inv(a), -e)
        r = self.one(); b = a
        while e:
            if e & 1: r = self.mul(r, b)
            b = self.mul(b, b); e >>= 1
        return r
    def taylor_F(self, pc, order):
        out = []
        for j in range(order + 1):
            acc = [GT()] * self.d
            for i in range(j, len(pc)):
                c = GT(pc[i])
                if c.is_zero(): continue
                tn = self.tpow(i - j); cc = c * comb(i, j)
                for m in range(self.d): acc[m] = acc[m] + cc * tn[m]
            out.append(tuple(acc))
        return out
    def taylor_K(self, pk, order):
        out = []
        for j in range(order + 1):
            acc = self.zero()
            for i in range(j, len(pk)):
                if self.is_zero(pk[i]): continue
                acc = self.add(acc, self.mul(pk[i], self.scal(comb(i, j), self.tpow(i - j))))
            out.append(acc)
        return out
    def ser_mul(self, A, B, order):
        out = []
        for n in range(order + 1):
            acc = self.zero()
            for i in range(n + 1):
                if i < len(A) and n - i < len(B) and not self.is_zero(A[i]) and not self.is_zero(B[n - i]): acc = self.add(acc, self.mul(A[i], B[n - i]))
            out.append(acc)
        return out
    def ser_inv(self, A, order):
        i0 = self.inv(A[0]); out = [i0]
        for n in range(1, order + 1):
            acc = self.zero()
            for i in range(1, n + 1):
                if i < len(A): acc = self.add(acc, self.mul(A[i], out[n - i]))
            out.append(self.scal(-1, self.mul(acc, i0)))
        return out
    def ser_pow(self, A, e, order):
        R = [self.one()] + [self.zero()] * order
        for _ in range(e): R = self.ser_mul(R, A, order)
        return R

def laurent_all(Ncoef, factors, lead, k):
    """R = N(y)/(lead * prod_i f_i^{m_i}); at a root rho of f_k (multiplicity m): returns (K, [h_0..h_{m-1}]) with
    R = sum_{j<m} h_j (y-rho)^{j-m} + regular, i.e. coefficient of (y-rho)^{-s} is h_{m-s}.  Elements of K = F0[y]/(f_k)."""
    fk, m = factors[k]; K = KF(fk); order = m - 1; d = K.d
    b = [None] * d
    b[d - 1] = K.const(K.f[d])
    for i in range(d - 1, 0, -1): b[i - 1] = K.add(K.const(K.f[i]), K.mul(K.tpow(1), b[i]))
    # f_k(y) = (y - rho) * B(y),  B(y) = sum_i b_i (y-rho)^i ... b above is B's coefficients in powers of y at... (convention: taylor_K(b) below)
    num = K.taylor_F(Ncoef, order)
    den = [K.const(lead)] + [K.zero()] * order
    den = K.ser_mul(den, K.ser_pow(K.taylor_K(b, order), m, order), order)
    for i, (fi, mi) in enumerate(factors):
        if i == k: continue
        den = K.ser_mul(den, K.ser_pow(K.taylor_F(fi, order), mi, order), order)
    h = K.ser_mul(num, K.ser_inv(den, order), order)
    return K, h

# ---------------------------------------------------------------- points and the tensor ring
class Points:
    """registry of fibration points.  A point id is ('Q', GT) or ('A', letter_idx, slot).  ZERO = ('Q', GT(0))."""
    def __init__(self, letter_xpolys):
        """letter_xpolys: {letter_idx: [GT coefficients of x^k]} for every x-dependent alphabet letter"""
        self.L = {i: ptrim(c) for i, c in letter_xpolys.items()}
        self.deg = {i: len(c) - 1 for i, c in self.L.items()}
        self._K = {}; self._rules = {}
    ZERO = ("Q", GT())
    def points_of_letter(self, i):
        d = self.deg[i]
        if d == 1:
            a1, a0 = self.L[i][1], self.L[i][0]; return [("Q", -(a0 / a1))]
        return [("A", i, s) for s in range(d)]
    def K(self, i):
        if i not in self._K: self._K[i] = KF(self.L[i])
        return self._K[i]
    def is_zero_pt(self, p): return p[0] == "Q" and p[1].is_zero()

class TRing:
    """elements: dict {exps: GT} where exps is a tuple of (gen, e) pairs sorted, gen = (letter, slot); with slot reduction rules."""
    def __init__(self, pts: Points):
        self.P = pts; self._rule = {}; self._pw = {}
    # ---- basic constructors
    @staticmethod
    def const(c):
        c = GT(c); return {} if c.is_zero() else {(): c}
    @staticmethod
    def is_zero(e): return not e
    def of_point(self, p):
        if p[0] == "Q": return self.const(p[1])
        return self.reduce({(((p[1], p[2]), 1),): GT(1)})
    # ---- arithmetic
    @staticmethod
    def add(a, b, scale=None):
        out = dict(a)
        for k, v in b.items():
            if scale is not None: v = v * scale
            w = out.get(k)
            if w is None: out[k] = v
            else:
                w = w + v
                if w.is_zero(): del out[k]
                else: out[k] = w
        return out
    @staticmethod
    def scal(a, c):
        c = GT(c)
        if c.is_zero(): return {}
        return {k: v * c for k, v in a.items()}
    @staticmethod
    def neg(a): return {k: -v for k, v in a.items()}
    @staticmethod
    def _mulmono(k1, k2):
        d = dict(k1)
        for g, e in k2: d[g] = d.get(g, 0) + e
        return tuple(sorted((g, e) for g, e in d.items() if e))
    def mul(self, a, b):
        if not a or not b: return {}
        out = {}
        for k1, v1 in a.items():
            for k2, v2 in b.items():
                k = self._mulmono(k1, k2); v = v1 * v2; w = out.get(k)
                if w is None: out[k] = v
                else:
                    w = w + v
                    if w.is_zero(): del out[k]
                    else: out[k] = w
        return self.reduce(out)
    def pow(self, a, e):
        assert e >= 0
        r = self.const(1); b = a
        while e:
            if e & 1: r = self.mul(r, b)
            e >>= 1
            if e: b = self.mul(b, b)
        return r
    # ---- slot rules: for letter i, slot s: monic polynomial P_s(y_s) of degree deg-s with coefficients (TRing elements in slots < s)
    def rule(self, i, s):
        key = (i, s)
        if key in self._rule: return self._rule[key]
        d = self.P.deg[i]
        if s == 0:
            lc = self.P.L[i][d]; coeffs = [self.const(c / lc) for c in self.P.L[i]]        # monic
        else:
            prev = self.rule(i, s - 1)          # coefficients (TRing) of P_{s-1}(y), degree d-s+1, monic
            g = ((i, s - 1), 1); gy = {(g,): GT(1)}
            # synthetic division of P_{s-1}(y) by (y - y_{s-1}):  q_{m-1} = p_m ; q_{k-1} = p_k + y_{s-1} q_k
            m = len(prev) - 1; q = [None] * m; q[m - 1] = prev[m]
            for k in range(m - 1, 0, -1): q[k - 1] = self.add(prev[k], self.mul(gy, q[k]))
            coeffs = q
        self._rule[key] = coeffs
        return coeffs
    def reduce(self, a):
        """apply slot rules until every generator exponent is below its rule degree"""
        out = {}; work = list(a.items())
        while work:
            k, v = work.pop()
            bad = None
            for g, e in k:
                dg = len(self.rule(g[0], g[1])) - 1
                if e >= dg: bad = (g, e, dg); break
            if bad is None:
                w = out.get(k)
                if w is None: out[k] = v
                else:
                    w = w + v
                    if w.is_zero(): del out[k]
                    else: out[k] = w
                continue
            g, e, dg = bad
            # y_g^e = y_g^(e-dg) * y_g^dg,  y_g^dg = - sum_{k<dg} r_k y_g^k
            rest = tuple((gg, ee) for gg, ee in k if gg != g); base = {rest: v}
            if e - dg: base = {self._mulmono(rest, ((g, e - dg),)): v}
            R = self.rule(g[0], g[1])
            for kk in range(dg):
                rk = R[kk]
                if not rk: continue
                term = {}
                for k1, v1 in base.items():
                    mono = self._mulmono(k1, ((g, kk),)) if kk else k1
                    for k2, v2 in rk.items():
                        kk2 = self._mulmono(mono, k2); term[kk2] = term.get(kk2, GT()) + (v1 * v2 * (-1))
                for kk2, vv in term.items():
                    if not vv.is_zero(): work.append((kk2, vv))
        return out
    # ---- K-polynomials (elements of KF as coefficient tuples in y) evaluated at a generator
    def from_K(self, coeffs, i, s):
        g = ((i, s), 1); out = {}
        for n, c in enumerate(coeffs):
            c = GT(c)
            if c.is_zero(): continue
            out[(((i, s), n),) if n else ()] = c
        return self.reduce(out)
    # ---- inverses
    def inv_point(self, p, k=1):
        """1/p^k"""
        if p[0] == "Q": return self.const(p[1].inv() ** k)
        i, s = p[1], p[2]; K = self.P.K(i)
        yinv = K.pow(K.inv(K.tpow(1)), k)
        return self.from_K(yinv, i, s)
    def point_pow(self, p, k):
        if k == 0: return self.const(1)
        if k < 0: return self.inv_point(p, -k)
        if p[0] == "Q": return self.const(p[1] ** k)
        return self.pow(self.of_point(p), k)
    def diff(self, a, p):
        """a - p"""
        return self.add(self.of_point(a), self.neg(self.of_point(p)))
    def inv_diff(self, a, p, k=1):
        """1/(a - p)^k for distinct points a, p"""
        key = ("invdiff", a if a[0] == "A" else ("Q", hash(a[1]), str(a[1])), p if p[0] == "A" else ("Q", hash(p[1]), str(p[1])), k)
        if key in self._pw: return self._pw[key]
        if k > 1:
            one = self.inv_diff(a, p, 1); res = self.pow(one, k); self._pw[key] = res; return res
        if a[0] == "Q" and p[0] == "Q":
            d = a[1] - p[1]
            if d.is_zero(): raise ZeroDivisionError("inv_diff of equal rational points")
            res = self.const(d.inv())
        elif a[0] == "A" and p[0] == "Q":
            K = self.P.K(a[1]); el = K.sub(K.tpow(1), K.const(p[1])); res = self.from_K(K.inv(el), a[1], a[2])
        elif a[0] == "Q" and p[0] == "A":
            K = self.P.K(p[1]); el = K.sub(K.const(a[1]), K.tpow(1)); res = self.from_K(K.inv(el), p[1], p[2])
        else:
            ia, sa = a[1], a[2]; ip, sp_ = p[1], p[2]
            if ia != ip:
                # 1/(y_a - y_p) = H(y_a, y_p) / M(y_a),  M = min poly of p (letter ip),  H = (M(y_a) - M(y_p))/(y_a - y_p)
                M = self.P.L[ip]; Ka = self.P.K(ia)
                My = Ka.from_ypoly(M)                      # M(y_a) in K_a
                invMy = self.from_K(Ka.inv(My), ia, sa)
                H = {}
                ya = self.of_point(a); yp = self.of_point(p)
                pa = [self.const(1)]; pp = [self.const(1)]
                dM = len(M) - 1
                for n in range(1, dM): pa.append(self.mul(pa[-1], ya)); pp.append(self.mul(pp[-1], yp))
                for kdeg in range(1, dM + 1):
                    mk = M[kdeg]
                    if mk.is_zero(): continue
                    for i_ in range(kdeg):          # y_a^i y_p^(k-1-i)
                        term = self.mul(pa[i_], pp[kdeg - 1 - i_]); H = self.add(H, term, scale=mk)
                res = self.mul(H, invMy)
            else:
                assert sa != sp_, "inv_diff of the identical algebraic point"
                d = self.P.deg[ia]; L = self.P.L[ia]; lc = L[d]; K = self.P.K(ia)
                # f'(y_a) = lc * prod_{j != a}(y_a - y_j)  ->  1/(y_a - y_p) = lc * prod_{j != a, p}(y_a - y_j) / f'(y_a)
                fprime = [L[n] * n for n in range(1, d + 1)]
                inv_fp = self.from_K(K.inv(K.from_ypoly(fprime)), ia, sa)
                others = [s_ for s_ in range(d) if s_ not in (sa, sp_)]
                prod = self.const(lc)
                for s_ in others: prod = self.mul(prod, self.diff(a, ("A", ia, s_)))
                res = self.mul(prod, inv_fp)
        self._pw[key] = res
        return res
    # ---- (de)serialization and numeric evaluation
    @staticmethod
    def ser(e): return [[[[g[0], g[1], ex] for (g, ex) in k], v.ser()] for k, v in e.items()]
    @staticmethod
    def deser(s): return {tuple(((g0, g1), ex) for g0, g1, ex in k): GT.deser(v) for k, v in s}
    @staticmethod
    def evalf(e, genval, tval, mp):
        """numeric value: genval {(letter, slot): mpc root}, t rational"""
        tot = mp.mpc(0)
        for k, v in e.items():
            x = v.mp(tval, mp)
            for g, ex in k: x = x * genval[g] ** ex
            tot += x
        return tot

def selftest():
    import mpmath as mp
    mp.mp.dps = 40
    t = GT.t(); i = GT.i()
    a = (t * 3 + i * 2 + 1) / (t * t + 1); b = a.inv(); assert (a * b) == GT(1)
    # KF over F0: f = y^3 - (t+i) y - 2
    f = [GT(-2), -(t + i), GT(0), GT(1)]; K = KF(f)
    y = K.tpow(1); el = K.add(K.mul(y, y), K.const(t)); inv = K.inv(el); one = K.mul(el, inv); assert one == K.one(), one
    # numeric check of TRing.inv_diff between roots of two letters and between two roots of the same cubic
    P = Points({7: f, 9: [GT(3), t, GT(1)], 11: [t + 1, GT(2)]})      # cubic, quadratic y^2 + t y + 3, linear 2y + (t+1)
    R = TRing(P); tv = Fr(2, 7)
    def roots(i_):
        cs = [c.mp(tv, mp) for c in P.L[i_]]; return mp.polyroots(list(reversed(cs)), maxsteps=200, extraprec=200)
    gv = {}
    for i_ in (7, 9):
        for s, r in enumerate(roots(i_)): gv[(i_, s)] = r
    out = {}
    A = ("A", 7, 0); Bq = ("A", 9, 1); A2 = ("A", 7, 2); Q = P.points_of_letter(11)[0]
    for nm, (u, v) in {"cubic-quad": (A, Bq), "cubic-cubic(slots 0,2)": (A, A2), "cubic-rational": (A, Q), "rational-quad": (Q, Bq), "quad slots": (("A", 9, 0), Bq)}.items():
        e = R.inv_diff(u, v, 2); val = TRing.evalf(e, gv, tv, mp)
        uu = gv[(u[1], u[2])] if u[0] == "A" else u[1].mp(tv, mp); vv = gv[(v[1], v[2])] if v[0] == "A" else v[1].mp(tv, mp)
        out[nm] = float(abs(val - 1 / (uu - vv) ** 2))
    e = R.inv_point(A, 3); out["inv_point"] = float(abs(TRing.evalf(e, gv, tv, mp) - 1 / gv[(7, 0)] ** 3))
    # reduce: y0^5 y1^4 numerically
    e = R.reduce({(((7, 0), 5), ((7, 1), 4)): GT(1)}); out["reduce"] = float(abs(TRing.evalf(e, gv, tv, mp) - gv[(7, 0)] ** 5 * gv[(7, 1)] ** 4))
    # laurent_all vs numeric: R = (y^2 + t)/( (2y+t+1)^2 * f(y) ) at roots of f and of the linear factor
    N = [t, GT(0), GT(1)]; facs = [([t + 1, GT(2)], 2), (f, 1)]
    Kf, h = laurent_all(N, facs, GT(1), 1); res = TRing.evalf(R.from_K(h[0], 7, 0), gv, tv, mp)
    yv = gv[(7, 0)]; fpr = sum(c.mp(tv, mp) * n * yv ** (n - 1) for n, c in enumerate(f) if n); ref = (yv ** 2 + t.mp(tv, mp)) / ((2 * yv + (t + 1).mp(tv, mp)) ** 2 * fpr)
    out["laurent_residue_cubic"] = float(abs(res - ref))
    Kl, h2 = laurent_all(N, facs, GT(1), 0)      # double pole at the rational root: h2[0] -> coefficient of (y-r)^-2, h2[1] -> ^-1
    r = -((t + 1) / 2); rv = r.mp(tv, mp); fv = lambda yy: sum(c.mp(tv, mp) * yy ** n for n, c in enumerate(f))
    g = lambda yy: (yy ** 2 + t.mp(tv, mp)) / (4 * fv(yy))
    c2 = g(rv); c1 = mp.diff(g, rv)
    out["laurent_double_pole_c2"] = float(abs(h2[0][0].mp(tv, mp) - c2)); out["laurent_double_pole_c1"] = float(abs(h2[1][0].mp(tv, mp) - c1))
    return out
if __name__ == "__main__":
    import json; r = selftest(); print(json.dumps(r, indent=1)); print("PASS", all(v < 1e-25 for v in r.values()))

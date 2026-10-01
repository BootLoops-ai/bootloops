"""zip_num — ZIP_reg (step-1 regulator) numeric evaluator.

Semantics: coef x prod ZIP_reg(word) with leading-1 shuffle regularization
and Hoelder convolution at 1/2, generalized to arbitrary NEGATIVE rational
letters + zeros, with trailing-zero handling at argument 1.  This is the
TRUSTED step-1 semantics of the HyperFLINT `op=hyperflint` integrator
(verified correct in all 869 calls of the production pilot).  Convention
pinned by live probe of the binary (100 digits vs independent quadrature).

INVARIANT (do not regress): a FIXED series length is a
dps-independent accuracy ceiling when mapped Hoelder letters approach 1/2
— a silent-truncation bug.  G_half therefore computes
its series length PER WORD from the actual minimum letter magnitude; the
constructor deliberately accepts NO fixed-length override.

Paper refs (SubTropica paper, reference/ vendored copy): hyperlogarithm
definition eq (3.40); shuffle product eq (3.43).
"""
from fractions import Fraction as Fr
import mpmath as mp


def shuffle(u, v):
    """All interweavings of tuples u, v with multiplicity (paper eq 3.43)."""
    if not u:
        return {v: 1}
    if not v:
        return {u: 1}
    out = {}
    for w, c in shuffle(u[1:], v).items():
        out[(u[0],) + w] = out.get((u[0],) + w, 0) + c
    for w, c in shuffle(u, v[1:]).items():
        out[(v[0],) + w] = out.get((v[0],) + w, 0) + c
    return out


class ZipEval:
    """Evaluator for ZIP_reg(word) at the working mpmath precision.

    Series lengths are chosen adaptively per word (see module docstring);
    caches (REG/TM/G1M/ZM) are per-instance, so build one instance per
    working precision.
    """

    def __init__(self):
        self.REG = {}   # leading-1 shuffle-reg decompositions
        self.TM = {}    # G(...;1/2) values
        self.G1M = {}   # G(...;1) via Hoelder convolution
        self.ZM = {}    # ZIP_reg values

    def reg_leading_ones(self, w):
        """Shuffle-regularize leading 1-letters (drop the divergent
        H(1;1,...) pieces): standard shuffle regularization at the
        endpoint, cf. paper eq (3.43) usage."""
        if w in self.REG:
            return self.REG[w]
        one = Fr(1)
        k = 0
        while k < len(w) and w[k] == one:
            k += 1
        if k == 0:
            self.REG[w] = {w: Fr(1)}
            return self.REG[w]
        v = w[k:]
        if not v:
            self.REG[w] = {}
            return self.REG[w]
        sh = shuffle((one,) * k, v)
        acc = {}
        for u, c in sh.items():
            if u == w:
                continue
            for uu, cc in self.reg_leading_ones(u).items():
                acc[uu] = acc.get(uu, Fr(0)) - Fr(c) * cc
        self.REG[w] = {u: c for u, c in acc.items() if c}
        return self.REG[w]

    def strip_trailing_zeros(self, w):
        """G(w;1): trailing zeros -> log(1)=0, so only the zero-log part
        of the shuffle-reg decomposition survives."""
        if not w or w[-1] != 0:
            return {w: Fr(1)}
        if all(a == 0 for a in w):
            return {}
        u = w[:-1]
        ins = {}
        for i in range(len(u) + 1):
            cand = u[:i] + (Fr(0),) + u[i:]
            ins[cand] = ins.get(cand, 0) + 1
        t = ins.pop(w)
        out = {}
        # log(1)*G(u;1) term vanishes
        for cand, m in ins.items():
            for v, c in self.strip_trailing_zeros(cand).items():
                out[v] = out.get(v, Fr(0)) - Fr(m) * c / t
        return out

    def G_half(self, w):
        """G(w; 1/2) by Taylor recursion.  ADAPTIVE series length: keyed
        to the minimum letter magnitude of THIS word (see
        module docstring)."""
        if not w:
            return mp.mpf(1)
        if w in self.TM:
            return self.TM[w]
        assert w[-1] != 0, w
        amin = min(abs(a) for a in w if a != 0)
        assert amin > Fr(1, 2), ('G_half letter too small', w)
        ratio = Fr(1, 2) / amin
        NT = int((mp.mp.dps + 12) * mp.log(10) /
                 (-mp.log(mp.mpf(ratio.numerator) / ratio.denominator))) + 10
        b = mp.mpf(w[-1].numerator) / w[-1].denominator
        c = [mp.mpf(0)] * (NT + 1)
        ib = 1 / b
        p = mp.mpf(1)
        for n in range(1, NT + 1):
            p *= ib
            c[n] = -p / n
        for a in reversed(w[:-1]):
            if a == 0:
                c = [mp.mpf(0)] + [c[n] / n for n in range(1, NT + 1)]
            else:
                assert abs(a) > Fr(1, 2), ('inner letter too small', a)
                av = mp.mpf(a.numerator) / a.denominator
                ia = 1 / av
                cp = [mp.mpf(0)] * (NT + 1)
                s = mp.mpf(0)
                ap = mp.mpf(1)
                iap = mp.mpf(1)
                for j in range(1, NT + 1):
                    s += c[j - 1] * ap
                    iap *= ia
                    cp[j] = -s * iap / j
                    ap *= av
                c = cp
        x = mp.mpf(1) / 2
        val = mp.mpf(0)
        xp = mp.mpf(1)
        for n in range(1, NT + 1):
            xp *= x
            val += c[n] * xp
        self.TM[w] = val
        return val

    def G_one_conv(self, w):
        """G(w; 1) by the Hoelder convolution at 1/2 (path composition,
        cf. paper §3.3.2 / eq 3.43 machinery)."""
        if not w:
            return mp.mpf(1)
        if w in self.G1M:
            return self.G1M[w]
        assert w[0] != 1 and w[-1] != 0, w
        n = len(w)
        tot = mp.mpf(0)
        for j in range(n + 1):
            left = tuple(Fr(1) - a for a in reversed(w[:j]))
            right = w[j:]
            tot += (-1) ** j * self.G_half(left) * self.G_half(right)
        self.G1M[w] = tot
        return tot

    def G_one(self, w):
        tot = mp.mpf(0)
        for v, c0 in self.strip_trailing_zeros(tuple(w)).items():
            for u, c in self.reg_leading_ones(v).items():
                q = c0 * c
                tot += mp.mpf(q.numerator) / q.denominator * self.G_one_conv(u)
        return tot

    def zip_value(self, word):
        """ZIP_reg(word), word = iterable of rationals, each <= 0 (letters
        on the integration contour are REFUSED — on-contour semantics are
        exactly the class where upstream fibration_basis is unsound)."""
        word = tuple(Fr(x) for x in word)
        if word in self.ZM:
            return self.ZM[word]
        for a in word:
            assert a <= 0, ('on-contour ZIP letter', a)
        exps = [((), Fr(1))]
        for a in word:
            opts = [(Fr(1), Fr(-1))]
            if a != Fr(-1):
                opts.append((a / (1 + a), Fr(1)))
            exps = [(w + (l,), c * s) for (w, c) in exps for (l, s) in opts]
        tot = mp.mpf(0)
        for w, c in exps:
            tot += mp.mpf(c.numerator) / c.denominator * self.G_one(w)
        self.ZM[word] = tot
        return tot

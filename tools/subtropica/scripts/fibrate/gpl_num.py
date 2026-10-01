"""gpl_num — mpmath GPL evaluator for {0,1,2}-string-letter words at
rational argument.  Taylor recursion + V-W trailing-zero shuffle strip.

Controlled in the F^red1
pilot against the stored close10_wout Arb word tables (36 words < 1e-70
relative).  Used downstream of fibrate to evaluate assembled formulas;
inside fibrate.py it serves as an independent cross-check of the
generalized evaluator (GEvalX) in --selfcheck.

Paper ref: hyperlogarithm definition eq (3.40) with the trailing-zero
convention G(0^k; x) = log^k(x)/k!.
"""
from fractions import Fraction as Fr
import mpmath as mp


def strip_zeros(w):
    """G(w;x) = sum coeff * log(x)^k * G(v;x), v without trailing zeros
    (shuffle-regularized trailing-zero strip)."""
    if not w or w[-1] != '0':
        return {(0, w): Fr(1)}
    if all(c == '0' for c in w):
        k = len(w)
        f = Fr(1)
        for j in range(1, k + 1):
            f /= j
        return {(k, ()): f}
    u = w[:-1]
    ins = {}
    for i in range(len(u) + 1):
        cand = u[:i] + ('0',) + u[i:]
        ins[cand] = ins.get(cand, 0) + 1
    t = ins.pop(w)
    out = {}
    for (k, v), c in strip_zeros(u).items():
        out[(k + 1, v)] = out.get((k + 1, v), Fr(0)) + c / t
    for cand, m in ins.items():
        for (k, v), c in strip_zeros(cand).items():
            out[(k, v)] = out.get((k, v), Fr(0)) - Fr(m) * c / t
    return out


class GEval:
    """G(word; x) for words over string letters '0','1','2',... at a fixed
    rational argument 0 < x < min nonzero letter."""

    def __init__(self, x_fr):
        self.x = mp.mpf(x_fr.numerator) / x_fr.denominator
        self.logx = mp.log(self.x)
        self.N = int((mp.mp.dps + 12) * mp.log(10) / (-mp.log(self.x))) + 10
        self.memo = {}

    def series(self, w):
        if w in self.memo:
            return self.memo[w]
        N = self.N
        b = mp.mpf(int(w[-1]))
        c = [mp.mpf(0)] * (N + 1)
        ib = 1 / b
        p = mp.mpf(1)
        for n in range(1, N + 1):
            p *= ib
            c[n] = -p / n
        for a in reversed(w[:-1]):
            if a == '0':
                c = [mp.mpf(0)] + [c[n] / n for n in range(1, N + 1)]
            else:
                av = mp.mpf(int(a))
                ia = 1 / av
                cp = [mp.mpf(0)] * (N + 1)
                s = mp.mpf(0)
                apow = mp.mpf(1)
                iapow = mp.mpf(1)
                for j in range(1, N + 1):
                    s += c[j - 1] * apow
                    iapow *= ia
                    cp[j] = -s * iapow / j
                    apow *= av
                c = cp
        self.memo[w] = c
        return c

    def G_noz(self, w):
        if not w:
            return mp.mpf(1)
        c = self.series(w)
        val = mp.mpf(0)
        xp = mp.mpf(1)
        for n in range(1, self.N + 1):
            xp *= self.x
            val += c[n] * xp
        return val

    def G(self, w):
        tot = mp.mpf(0)
        for (k, v), q in strip_zeros(tuple(w)).items():
            tot += (mp.mpf(q.numerator) / q.denominator) * self.logx ** k \
                   * self.G_noz(v)
        return tot

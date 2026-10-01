#!/usr/bin/env python3
# ============================================================================
# TRUST vendored engine: bessel_oracle2.py (sha pinned in trust/_pins.py).
# ============================================================================
"""Oracle v2: sec39-class = 4-line equal-mass BANANA at P^2 = t = -1/3
(NOT q=0: G = e3(1+M e1) - t*e4, verified from build_G). Euclid p^2 = +1/3.

W_kira(nu; p) = (4pi)^{3d/2} (2pi)^{d/2} int_0^inf dr r^{d-1} ang(r) prod_i f_{nu_i}(r)
ang(r) = (pr)^{1-d/2} J_{d/2-1}(pr) = 2^{1-d/2} sum_j (-p^2/4)^j r^{2j}/(j! Gamma(d/2+j))
f_nu as in v1. I_kira = (-1)^w * Euclid. Tadpole^3 column = T(1)^3 (G checked).
Continuation: series head on [0,1/2] (exact Fraction exponents) + numeric tail.
Validation: (a) series-vs-direct at d=5/2; (b) dW/dM = -4 W(1,1,1,2) via
central difference (exact identity, any p).
"""
import json, sys
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 70
NTRUNC = 80
RCUT = mp.mpf(1)/2

def frac_mp(fr): return mp.mpf(fr.numerator)/mp.mpf(fr.denominator)

class Oracle:
    def __init__(self, d_frac, M, p2):
        self.dF = d_frac; self.d = frac_mp(d_frac)
        self.M = mp.mpf(M); self.p2 = mp.mpf(p2); self.p = mp.sqrt(self.p2)

    def tad(self, nu):
        return mp.gamma(nu - self.d/2)/mp.gamma(nu) * self.M**(self.d/2 - nu)

    def f_bessel(self, nu, r):
        a = nu - self.d/2
        return 2*(4*mp.pi)**(-self.d/2)/mp.gamma(nu) \
               * (r/(2*mp.sqrt(self.M)))**a * mp.besselk(a, mp.sqrt(self.M)*r)

    def ang(self, r):
        return (self.p*r)**(1-self.d/2) * mp.besselj(self.d/2-1, self.p*r)

    def f_series(self, nu):
        aF = Fraction(nu) - self.dF/2; a = frac_mp(aF)
        c = 2*(4*mp.pi)**(-self.d/2)/mp.gamma(nu) * mp.pi/(2*mp.sin(a*mp.pi))
        s0 = [c*self.M**(-a)/(mp.factorial(j)*mp.gamma(-a+j+1))*(self.M/4)**j
              for j in range(NTRUNC)]
        s1 = [-c*mp.mpf(4)**(-a)/(mp.factorial(j)*mp.gamma(a+j+1))*(self.M/4)**j
              for j in range(NTRUNC)]
        return [(Fraction(0), s0), (2*aF, s1)]

    def ang_series(self):
        c = mp.mpf(2)**(1-self.d/2)
        s = [c*(-self.p2/4)**j/(mp.factorial(j)*mp.gamma(self.d/2+j))
             for j in range(NTRUNC)]
        return [(Fraction(0), s)]

    def W(self, nus, direct=False):
        pref = (4*mp.pi)**(3*self.d/2) * (2*mp.pi)**(self.d/2)
        def igr(r):
            v = r**(self.d-1) * self.ang(r)
            for nu in nus: v *= self.f_bessel(nu, r)
            return v
        if direct:
            return pref*mp.quad(igr, [0, RCUT, 2, 8, 30, mp.inf])
        factors = [self.ang_series()] + [self.f_series(nu) for nu in nus]
        prod = {Fraction(0): [mp.mpf(1)] + [mp.mpf(0)]*(NTRUNC-1)}
        for fs in factors:
            new = {}
            for pe, pc in prod.items():
                for be, bc in fs:
                    e = pe + be
                    conv = new.setdefault(e, [mp.mpf(0)]*NTRUNC)
                    for i in range(NTRUNC):
                        if pc[i] == 0: continue
                        for j in range(NTRUNC - i):
                            if bc[j] == 0: continue
                            conv[i+j] += pc[i]*bc[j]
            prod = new
        head = mp.mpf(0); base = self.dF - 1
        for e, cs in prod.items():
            for j, cj in enumerate(cs):
                b = base + e + 2*j
                bm = frac_mp(b)
                head += cj * RCUT**(bm+1)/(bm+1)
        tail = mp.quad(igr, [RCUT, 2, 8, 30, mp.inf])
        return pref*(head + tail)

def main():
    d0F = Fraction(97,23); M = mp.mpf(12)/7; p2 = mp.mpf(1)/3
    # (a) strip selftest
    ot = Oracle(Fraction(5,2), M, p2)
    for nus in [(1,1,1,1), (1,1,1,2)]:
        a = ot.W(nus); b = ot.W(nus, direct=True)
        dig = float(-mp.log10(abs(a-b)/abs(b)))
        print(f"[selftest d=5/2, p2=1/3] W{nus}: {dig:.1f}d ({'PASS' if dig>30 else 'FAIL'})")
    o = Oracle(d0F, M, p2)
    # (b) dW/dM = -4 W(...,2) exact identity, central diff
    h = mp.mpf(10)**(-25)
    op = Oracle(d0F, M+h, p2); om = Oracle(d0F, M-h, p2)
    lhs = (op.W((1,1,1,1)) - om.W((1,1,1,1)))/(2*h)
    rhs = -4*o.W((1,1,1,2))
    dig = float(-mp.log10(abs(lhs-rhs)/abs(rhs)))
    print(f"[massderiv] dW/dM vs -4W(1,1,1,2): {dig:.1f}d ({'PASS' if dig>15 else 'FAIL'})")

    Wcache = {}
    def col_value(name):
        if name == '((1, 2, 3), (0, 0, 0))':
            return (-1)**3 * o.tad(1)**3
        m = eval(name.split('), ')[1].rstrip(')') + ')')
        nus = tuple(sorted(mi+1 for mi in m))
        if nus not in Wcache: Wcache[nus] = o.W(nus)
        return (-1)**sum(nus) * Wcache[nus]
    def jfac(name):
        m = (0,0,0) if name == '((1, 2, 3), (0, 0, 0))' else eval(name.split('), ')[1].rstrip(')') + ')')
        w = sum(m) + len(m)
        c = mp.mpf(1)
        for mi in m: c *= mp.factorial(mi)
        for j in range(1, w+1): c /= (2*o.d - j)
        return c * (-1)**w
    def dot(vec, space='I'):
        tot = mp.mpf(0); scale = mp.mpf(0)
        for name, c in vec.items():
            cf = Fraction(c)
            cv = mp.mpf(cf.numerator)/mp.mpf(cf.denominator)
            v = col_value(name)
            if space == 'J': v = jfac(name)*v
            tot += cv*v; scale += abs(cv*v)
        return tot, scale

    r = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'sec63_final.json'))
    print("\n[1] violation rows in J-space (must vanish if true):")
    seen = set()
    for i, vr in enumerate(r['viol_rows_named']):
        key = tuple(sorted(vr.items()))
        if key in seen: continue
        seen.add(key)
        tot, scale = dot(vr, space='J')
        dig = float(-mp.log10(abs(tot)/scale)) if tot != 0 else 99.0
        print(f"  viol[{i}]: rel.residual = 10^-{dig:.1f}  ({'PASS' if dig>25 else 'FAIL'})")
    print("\n[2] Delta = mine - kira, I-space (vanishes iff both valid):")
    for d in r['details_sigma_-1']:
        if d['match']: continue
        delta = {}
        for name, mine, kira in d['diffs']:
            delta[name] = str(Fraction(mine) - Fraction(kira))
        tot, scale = dot(delta)
        dig = float(-mp.log10(abs(tot)/scale)) if tot != 0 else 99.0
        print(f"  target {d['target'][:6]}: rel.residual = 10^-{dig:.1f}  ({'PASS' if dig>25 else 'FAIL'})")
    print("\n[3] kira side alone vs mine alone (which one is a valid reduction?):")
    print("  (residual of target reduction cannot be tested without target value;")
    print("   only Delta is testable — [2] is the decisive gate.)")

if __name__ == '__main__':
    main()

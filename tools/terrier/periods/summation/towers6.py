#!/usr/bin/env python3
"""towers6.py — B3: 6-var rho-jet Gamma-series tower evaluator for HV4 (5.18).
Objects (jet depth <= 4, 210 monomials):
RHO towers  T_j(phi) = sum_{|k|<=M} c(k) E_j(k) phi^k,   c(k)=multinom(|k|;k)^2,
  E(k;rho) = exp( 2 sum_{m=1..4} s_m [ H^(m)_{|k|} r^m - sum_I H^(m)_{k_I} rho_I^m ] ),
  s_m = (-1)^(m-1)/m, r = sum_I rho_I, truncated at total rho-degree 4.
  (Pochhammer: Gamma(1+n+x)/[Gamma(1+n)Gamma(1+x)] = exp(sum_m s_m H^(m)_n x^m); EXACT.
   These are the A(phi) towers of Pi = exp(t.N) A(phi): phi^rho prefactor lives in exp(t.N).)
THETA jets  a_alpha(phi) = theta^alpha Pi^0 = sum_{|k|<=M} c(k) prod_I k_I^alpha_I phi^k
  (the B7_CAUCHY.md interface: certified center balls for a_alpha, |alpha| <= 4).
EVALUATOR (separable GF form; EXACTLY equal to the lattice sum, no approximation):
  g_I(b) = phi_I^b/(b!)^2 * jet_I(b);  G_I(x) = sum_b g_I(b) x^b;
  P = prod_I G_I |_{x^M, deg<=4};  T = sum_n (n!)^2 F_n(r) [x^n] P  (F_n = 1 theta mode).
Rings: exact (Fraction / CQ complex-rational) and flint arb/acb balls (prec set by caller).
Tail law (entropy majorant, README.md in this directory): see tail_theta /
tail_rho below."""
from fractions import Fraction as Fr
from math import factorial

H, JDEG = 6, 4
PB = [5 ** i for i in range(H)]          # packed base-5 monomials (no carries)
_MONS = []                                # all monos deg<=4 as packed ints


def _gen_mons(pos, left, cur):
    if pos == H:
        _MONS.append(cur)
        return
    for d in range(left + 1):
        _gen_mons(pos + 1, left - d, cur + d * PB[pos])


_gen_mons(0, JDEG, 0)
DEG = {}
for m in _MONS:
    t, x = 0, m
    for i in range(H):
        t += x % 5
        x //= 5
    DEG[m] = t
NMON = len(_MONS)                         # 210


def unpack(m):
    return tuple((m // PB[i]) % 5 for i in range(H))


def jetmul(A, B):
    """product of jet dicts {packed_mono: coeff}, total degree capped at 4."""
    out = {}
    for ma, ca in A.items():
        da = DEG[ma]
        for mb, cb in B.items():
            if da + DEG[mb] <= JDEG:
                m = ma + mb
                v = ca * cb
                if m in out:
                    out[m] = out[m] + v
                else:
                    out[m] = v
    return out


def jet_exp_uni(coeffs):
    """exp of sum_{m=1..4} coeffs[m] t^m, truncated deg<=4 -> [a0..a4] exact."""
    a = [Fr(1), Fr(0), Fr(0), Fr(0), Fr(0)]
    pw = [Fr(1), Fr(0), Fr(0), Fr(0), Fr(0)]     # Q^j / j!
    for j in range(1, 5):
        nxt = [Fr(0)] * 5
        for d1 in range(5):
            if pw[d1] == 0:
                continue
            for m in range(1, 5):
                if d1 + m <= 4 and coeffs[m]:
                    nxt[d1 + m] += pw[d1] * coeffs[m]
        pw = [x / j for x in nxt]
        for d in range(5):
            a[d] += pw[d]
    return a


# harmonic tables H^(m)_n, m=1..4 (exact Fractions)
def harmonic_tables(N):
    Ht = {m: [Fr(0)] * (N + 1) for m in range(1, 5)}
    for n in range(1, N + 1):
        for m in range(1, 5):
            Ht[m][n] = Ht[m][n - 1] + Fr(1, n ** m)
    return Ht


class CQ:
    """exact complex rational (re, im Fractions); ring element for gf_eval."""
    __slots__ = ("re", "im")

    def __init__(self, re, im=Fr(0)):
        self.re, self.im = Fr(re), Fr(im)

    def __add__(self, o):
        return CQ(self.re + o.re, self.im + o.im)

    def __mul__(self, o):
        if isinstance(o, CQ):
            return CQ(self.re * o.re - self.im * o.im,
                      self.re * o.im + self.im * o.re)
        return CQ(self.re * o, self.im * o)

    def __eq__(self, o):
        return self.re == o.re and self.im == o.im

    __rmul__ = __mul__

    def __repr__(self):
        return f"CQ({self.re},{self.im})"


def s_m(m):
    return Fr((-1) ** (m - 1), m)


def axis_jet_coeffs(b, Ht, sign):
    """[a0..a4](Fr): exp(sign * 2 sum_m s_m H^(m)_b t^m) truncated deg 4."""
    return jet_exp_uni([None] + [sign * 2 * s_m(m) * Ht[m][b]
                                 for m in range(1, 5)])


def rpow_expansions():
    """r^d = (rho_1+..+rho_6)^d as jet dicts with integer multinomials, d=0..4."""
    out = [{0: 1}]
    r1 = {PB[i]: 1 for i in range(H)}
    for _ in range(4):
        out.append(jetmul(out[-1], r1))
    return out


def g_tables_exact(phis, M, mode, Ht):
    """per-axis exact CQ jet coefficients g[I][b] = {d*PB[I]: CQ}, b=0..M.
    mode 'rho': jet = exp(-2 sum s_m H^(m)_b rho_I^m); 'theta': jet_d = b^d."""
    g = []
    for I in range(H):
        col, pw = [], CQ(1)
        for b in range(M + 1):
            base = pw * Fr(1, factorial(b) ** 2)
            if mode == "rho":
                cf = axis_jet_coeffs(b, Ht, -1)
                col.append({d * PB[I]: base * cf[d]
                            for d in range(5) if cf[d]})
            else:
                col.append({d * PB[I]: base * (b ** d)
                            for d in range(5) if b ** d or d == 0})
            pw = pw * phis[I]
        g.append(col)
    return g


def contraction_weights_exact(M, mode, Ht):
    """w[n][d] = (n!)^2 * F_n coeff of r^d (mode rho); theta: w[n][0] only."""
    out = []
    for n in range(M + 1):
        w2 = Fr(factorial(n) ** 2)
        if mode == "rho":
            f = axis_jet_coeffs(n, Ht, +1)
            out.append([w2 * f[d] for d in range(5)])
        else:
            out.append([w2, 0, 0, 0, 0])
    return out


def gf_eval(gtab, wts, M, rpow):
    """separable GF pass. gtab/wts/rpow already ring-converted (wts entry None
    = exact zero, skip). Returns jet dict {packed_mono: ring}."""
    P = [dict(gtab[0][b]) for b in range(M + 1)]
    for I in range(1, H):
        col, newP = gtab[I], []
        for n in range(M + 1):
            acc = {}
            for a in range(n + 1):
                gb = col[n - a]
                for ma, ca in P[a].items():
                    da = DEG[ma]
                    for mb, cb in gb.items():
                        if da + DEG[mb] <= JDEG:
                            m = ma + mb
                            v = ca * cb
                            if m in acc:
                                acc[m] = acc[m] + v
                            else:
                                acc[m] = v
            newP.append(acc)
        P = newP
    S = [None] * 5
    for n in range(M + 1):
        for d in range(5):
            w = wts[n][d]
            if w is None:
                continue
            tgt = S[d] if S[d] is not None else {}
            for m, c in P[n].items():
                if m in tgt:
                    tgt[m] = tgt[m] + c * w
                else:
                    tgt[m] = c * w
            S[d] = tgt
    T = {}
    for d in range(5):
        if S[d] is None:
            continue
        for m, c in jetmul(rpow[d], S[d]).items():
            if m in T:
                T[m] = T[m] + c
            else:
                T[m] = c
    return T


def lattice_eval_exact(phis, M, mode, Ht):
    """independent direct-lattice exact path (gate engine): enumerate |k|<=M,
    E(k) = truncated jet-exp of the FULL 6-var Q_k (generic exp, no per-axis
    factorization); theta mode: plain monomial weights. Exact CQ output."""
    RP = rpow_expansions()
    T = {}
    stack = [((), 0)]
    while stack:
        pre, n = stack.pop()
        if len(pre) == H - 1:
            for kl in range(M - n + 1):
                k = pre + (kl,)
                nn = n + kl
                den = 1
                for x in k:
                    den *= factorial(x) ** 2
                w = CQ(Fr(factorial(nn) ** 2, den))
                for I in range(H):
                    for _ in range(k[I]):
                        w = w * phis[I]
                if mode == "theta":
                    for mo in _MONS:
                        u, wt = unpack(mo), 1
                        for I in range(H):
                            wt *= k[I] ** u[I]
                        if wt:
                            T[mo] = T.get(mo, CQ(0)) + w * wt
                    continue
                Q = {}
                for m in range(1, 5):
                    cm = 2 * s_m(m) * Ht[m][nn]
                    for mo, mult in RP[m].items():
                        Q[mo] = Q.get(mo, Fr(0)) + cm * mult
                    for I in range(H):
                        Q[m * PB[I]] = (Q.get(m * PB[I], Fr(0))
                                        - 2 * s_m(m) * Ht[m][k[I]])
                Q = {mo: c for mo, c in Q.items() if c}
                E, pw = {0: Fr(1)}, {0: Fr(1)}
                for j in range(1, 5):
                    pw = jetmul(pw, Q)
                    pw = {mo: c / j for mo, c in pw.items()}
                    for mo, c in pw.items():
                        E[mo] = E.get(mo, Fr(0)) + c
                for mo, c in E.items():
                    if c:
                        T[mo] = T.get(mo, CQ(0)) + w * c
        else:
            for kx in range(M - n + 1):
                stack.append((pre + (kx,), n + kx))
    return T


def J_rho_exact(m, Ht):
    """J_d(m), d=0..4: max |j|=d coeff of exp(Qbar_m), Qbar = 2 H_m *
    sum_{m'=1..4} [ (sum rho)^m' + sum_I rho_I^m' ] / m'  (majorant; EXACT)."""
    RP = rpow_expansions()
    Hm = Ht[1][m]
    Q = {}
    for mp in range(1, 5):
        c = 2 * Hm / mp
        for mo, mult in RP[mp].items():
            Q[mo] = Q.get(mo, Fr(0)) + c * mult
        for I in range(H):
            Q[mp * PB[I]] = Q.get(mp * PB[I], Fr(0)) + c
    E, pw = {0: Fr(1)}, {0: Fr(1)}
    for j in range(1, 5):
        pw = jetmul(pw, Q)
        pw = {mo: c / j for mo, c in pw.items()}
        for mo, c in pw.items():
            E[mo] = E.get(mo, Fr(0)) + c
    J = [Fr(0)] * 5
    for mo, c in E.items():
        d = DEG[mo]
        if c > J[d]:
            J[d] = c
    return J


def tails_at(s2, M, Ht, arb):
    """certified tail upper bounds per jet degree d=0..4, BOTH modes, at
    entropy scale s^2 = s2 (arb ball; use .upper() outside). Returns
    (tail_theta[d], tail_rho[d]) arb lists; None entry = hypothesis s2*e^eps<1
    FAILED (fail-closed: caller must treat as NO certificate).
    theta: sum_{m>M} m^d x^m <= (M+1)^d x^(M+1) / (1 - x e^(d/(M+1)))
    rho:   sum_{m>M} J_d(m) x^m <= J_d(M+1) x^(M+1) / (1 - x e^(d/(H_(M+1)(M+1))))
    [J_d poly deg<=d in H_m, nonneg coeffs => J_d(m) <= J_d(M+1)(H_m/H_{M+1})^d;
     H_m <= H_{M+1}+ln(m/(M+1)) and 1+x<=e^x => (H_m/H_{M+1})^d <= e^(d(m-M-1)/(H_{M+1}(M+1)))]"""
    J = J_rho_exact(M + 1, Ht)
    HM1 = Ht[1][M + 1]
    one = arb(1)
    xt = s2 ** (M + 1)
    th, rh = [], []
    for d in range(5):
        c = s2 * (arb(d) / (M + 1)).exp()
        th.append((arb(M + 1) ** d) * xt / (one - c) if c < 1 else None)
        c2 = s2 * (arb(d) / (arb(HM1.numerator) / HM1.denominator * (M + 1))).exp()
        rh.append((arb(J[d].numerator) / J[d].denominator) * xt / (one - c2)
                  if c2 < 1 else None)
    return th, rh


def conv_arb_maker(arb):
    def cv(x):
        if isinstance(x, CQ):
            assert x.im == 0, "arb ring needs real point"
            x = x.re
        return arb(x.numerator) / x.denominator
    return cv


def conv_acb_maker(arb, acb):
    def cv(x):
        if isinstance(x, CQ):
            return acb(arb(x.re.numerator) / x.re.denominator,
                       arb(x.im.numerator) / x.im.denominator)
        return acb(arb(x.numerator) / x.denominator)
    return cv


def ring_tables(gtab, wts, cv):
    """convert exact g tables + contraction weights through cv; 0 -> None."""
    g2 = [[{m: cv(c) for m, c in cell.items()} for cell in col]
          for col in gtab]
    w2 = [[(cv(w) if w else None) for w in row] for row in wts]
    return g2, w2


def sq_modulus(phi):
    """|phi|^2 exact Fraction of a CQ."""
    return phi.re * phi.re + phi.im * phi.im

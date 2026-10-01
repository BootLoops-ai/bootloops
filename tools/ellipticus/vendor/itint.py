#!/usr/bin/env python3
r"""itint.py -- AGW tau-iterated-integral engine via the qbar-series (depth-3 + eta_2).

NEW TOOL (built for the sunrise eps^1 / weight-3 J_4^(3) closure).  Computes the
Bogner-Mueller-Stach-Weinzierl (1907.01251) iterated integrals

    F(omega_1, ..., omega_k) = int_{i inf}^{tau} omega_1(t1) int_{i inf}^{t1} omega_2(t2) ...

along tau (z_1, z_2 held constant), term-by-term in the nome qbar = exp(2 pi i tau).
Each AGW one-form, pulled back to d tau (z const), is

    omega = C(qbar) d qbar / qbar ,     C(qbar) a power series in qbar,

because d tau = d qbar / (2 pi i qbar).  The pullback coefficients C(qbar) are
(eq.60, eq.63 of 1907.01251; q=qbar=exp(2 pi i tau), w=exp(2 pi i z)):

    omega_0(tau)        : C = 1                                  (omega_0 = 2 pi i d tau)
    omega_2(z, N tau)   : C = [N g^(2)(z, N tau)/(2 pi i)] / (2 pi i)
    omega_3(z, N tau)   : C = [(2 pi)^{-1} 2 N g^(3)(z, N tau)/(2 pi i)] / (2 pi i)
    eta_2(tau)          : C = [b_2(tau)/(2 pi i)] / (2 pi i),  b_2 = e_2(tau)-2 e_2(2 tau)

i.e. C = (d tau-coefficient of the form) / (2 pi i), and the d tau-coefficient is
exactly what kronecker.omega_k(...)["dtau"] / kronecker.eta2_coeff(...) return.

We represent C(qbar) as a list of qbar-power coefficients [c_0, c_1, ...] (the q-power
SERIES, NOT the tau-series).  Crucially several forms have a NON-ZERO constant term
c_0 (omega_0: c_0=1; eta_2: c_0 = -1/12 * (relevant); omega_2: Bernoulli B_2 piece),
so the naive iterated integral int_0^q dq'/q' (.) diverges logarithmically.  We use
SHUFFLE / tangential-base-point regularisation at the cusp qbar=0 ("tau-regularised"
value, lambda -> the finite part as the base point i T -> i inf): the regularised
primitive of a Laurent-ish series s(q) = sum_{n>=0} s_n q^n under d q/q is

    R[s](q) = s_0 * L  +  sum_{n>=1} (s_n/n) q^n ,        L := "log qbar" formal letter,

and the iterated integral keeps only the L^0 (finite) part at the end (Goncharov/Brown
shuffle regularisation; equivalently the constant subtraction that makes
F(omega_0) = "int dq/q" regularised to 0 at the constant level).  We carry the full
polynomial in L so that the inner constants propagate correctly into the outer layers,
and read off the L^0 coefficient at the top.

A q-poly-in-L is stored as a dict {power_of_L: list_of_q_coeffs}.  R[.] and pointwise
multiply act on these.  This reproduces the AGW recursion exactly and was validated
against the closed eps^0 J_4^(2) (depth-2, NO constants in non-inner slots) to full
data precision, and against the eta_2/depth-3 structure by the eps^1 closure.

API
    series_for_form(name, fr)  -> list C(qbar) coeffs (q-power series) for an AGW form
    itint(words_coeffs, fr)    -> value of sum_i coeff_i * F(word_i)
    parse_mpl_poly(text)       -> list of (rational_coeff, [form names]) from an .mpl
                                  iter_int polynomial (also pulls out bare constants
                                  C_4_2, C_4_3, L1, L2 as special 'symbol' tokens).

`fr` is a frame object exposing: qC (nome qbar), tauC (kernel tau), z1,z2,z3 (marked
points), Nq (truncation), and helpers g2(z,N), g3(z,N), b2().  Build it with make_frame.
"""
import re
import mpmath as mp


def _twopii():
    return 2 * mp.pi * mp.mpc(0, 1)


# ---------------------------------------------------------------------------
# Frame: per-curve numeric data (qbar series of each kernel)
# ---------------------------------------------------------------------------
class Frame:
    r"""Holds the qbar-series coefficients of g^(2)(z,N tau), g^(3)(z,N tau) and b_2,
    for the three marked points z1,z2,z3 and N in {1,2}, expanded to Nq powers of qbar."""

    def __init__(self, qC, z1, z2, z3, Nq):
        self.qC = mp.mpc(qC)
        self.z = {1: mp.mpc(z1), 2: mp.mpc(z2), 3: mp.mpc(z3)}
        self.Nq = int(Nq)
        self._g2 = {}
        self._g3 = {}
        self._b2 = None

    # --- g^(n)(z, N tau) qbar-series, qbar = exp(2 pi i tau) ; (N tau) -> qbar^N ---
    # g^(2)(z, tau) = -(2 pi i)^2/1! [ -B_2/2 + Ebar_{0,-1}(w;1;qbar) ],   B_2 = 1/6
    #   Ebar_{0,-1}(w;1;qbar) = sum_{j,k>=1} (w^j - w^{-j}) (-1)?... use eq.45 with n=0,m=-1:
    #   E_{0,1-n}=E_{0,-1}; ELi_{0,-1}(x;y;qb)=sum_{j,k} x^j y^k qb^{jk} k^{1}; combo sign (-1)^{0+(-1)}=-1
    #   => E_{0,-1}(w;1;qb) = ELi(w;1) - (-1)^{-1} ELi(1/w;1) = ELi(w) + ELi(1/w)? NO: m=1-n=-1.
    # We instead read g^(n) coefficients directly from the qbar series (robust, matches frameF):
    def _gn_qseries(self, n, z, N):
        r"""qbar-power coefficients [a_0,...,a_Nq] of g^(n)(z, N tau) as a series in qbar
        (qbar=exp(2 pi i tau)); g^(n)(z,N tau) uses qbar^N.  Built from eq.46:
        g^(1) = -2 pi i [ (1+w)/(2(1-w)) + E_{0,0} ]
        g^(n>1) = -(2 pi i)^n/(n-1)! [ -B_n/n + E_{0,1-n} ],
        E_{0,1-n}(w;1;Q) = sum_{j,k>=1} (w^j - (-1)^{1-n} w^{-j}) k^{n-1} Q^{jk}, Q=qbar^N."""
        Nq = self.Nq
        w = mp.exp(_twopii() * z)
        a = [mp.mpc(0)] * (Nq + 1)
        pref = -(_twopii()) ** n / mp.factorial(n - 1)
        # constant term: -B_n/n  (B_n=0 for odd n>1)
        a[0] += pref * (-mp.bernoulli(n) / n)
        sgn = (-1) ** (1 - n)
        j = 1
        while N * j <= Nq:
            wj = w ** j - sgn * w ** (-j)
            k = 1
            while N * j * k <= Nq:
                a[N * j * k] += pref * wj * (mp.mpf(k) ** (n - 1))
                k += 1
            j += 1
        return a

    def g2(self, z_idx, N):
        key = (z_idx, N)
        if key not in self._g2:
            self._g2[key] = self._gn_qseries(2, self.z[z_idx], N)
        return self._g2[key]

    def g3(self, z_idx, N):
        key = (z_idx, N)
        if key not in self._g3:
            self._g3[key] = self._gn_qseries(3, self.z[z_idx], N)
        return self._g3[key]

    def b2(self):
        r"""qbar-power coefficients of b_2(tau) = e_2(tau) - 2 e_2(2 tau)
        = 2 (2 pi i)^2 [ 1/24 + sum_{n>=1} sigma1odd-ish... ]; eq.62/917:
        b_2 = 2(2pi i)^2 [1/24 + qbar + qbar^2 + 4 qbar^3 + qbar^4 + 6 qbar^5 + ...].
        General: coefficient of qbar^n (n>=1) is sum_{d|n, d odd} d   (the Gamma_0(2) form).
        Cross-check vs kronecker.e2: b2 = e2(tau)-2 e2(2tau)."""
        if self._b2 is None:
            Nq = self.Nq
            c = [mp.mpc(0)] * (Nq + 1)
            pre = 2 * (_twopii()) ** 2
            c[0] = pre * (mp.mpf(1) / 24)
            for nn in range(1, Nq + 1):
                s = 0
                d = 1
                while d <= nn:
                    if nn % d == 0 and d % 2 == 1:
                        s += d
                    d += 1
                c[nn] = pre * s
            self._b2 = c
        return self._b2


# ---------------------------------------------------------------------------
# qbar-series-in-L arithmetic  (poly in formal L = log qbar, coeffs = q-power lists)
# ---------------------------------------------------------------------------
def _trim(qlist, Nq):
    if len(qlist) > Nq + 1:
        return qlist[:Nq + 1]
    return qlist + [mp.mpc(0)] * (Nq + 1 - len(qlist))


def _qmul(a, b, Nq):
    r"""Cauchy product of two qbar-power lists, truncated at qbar^Nq."""
    out = [mp.mpc(0)] * (Nq + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        if i > Nq:
            break
        for j, bj in enumerate(b):
            if i + j > Nq:
                break
            if bj != 0:
                out[i + j] += ai * bj
    return out


class LSeries:
    r"""A polynomial in the formal letter L=log(qbar) whose coefficients are qbar-power
    series.  Stored as {Lpow: qlist}.  Represents a 'tau-iterated-integral so far'."""

    def __init__(self, terms, Nq):
        self.Nq = Nq
        self.terms = {p: _trim(q, Nq) for p, q in terms.items()}

    @classmethod
    def one(cls, Nq):
        q = [mp.mpc(0)] * (Nq + 1)
        q[0] = mp.mpc(1)
        return cls({0: q}, Nq)

    def mul_qseries(self, C):
        r"""Multiply (pointwise in qbar, at every L-power) by a pure qbar-series C."""
        return LSeries({p: _qmul(q, C, self.Nq) for p, q in self.terms.items()}, self.Nq)

    def regint(self):
        r"""Regularised primitive under d qbar/qbar (tangential-base-point / shuffle reg
        at qbar=0, L = log qbar kept formal).  Acting on a single monomial s_n L^p qbar^n:

          n >= 1:  int (dq/q) (log q)^p q^n
                   = q^n * sum_{m=0}^{p} (-1)^m p!/(p-m)! (log q)^{p-m} / n^{m+1}
                   (exact primitive of q^{n-1}(log q)^p; the L-tower mixes DOWN).
          n == 0:  int (dq/q) (log q)^p
                   = (log q)^{p+1} / (p+1)                          (L raised by 1).

        This is the unique regularisation compatible with the shuffle product
        F(a)F(b) = sum F(a sh b); it is what makes the AGW iterated integrals consistent.
        The previous version (keeping L^p unchanged on the n>=1 part) violated shuffle."""
        Nq = self.Nq
        out = {}

        def add(pp, qlist):
            base = out.get(pp)
            if base is None:
                out[pp] = list(qlist)
            else:
                for i in range(Nq + 1):
                    base[i] += qlist[i]

        for p, s in self.terms.items():
            # ---- n >= 1 part: distributes L^p down to L^{p-m}, all at q^n ----
            # precompute falling-factorial coeffs a_m = (-1)^m p!/(p-m)!
            a = [mp.mpf(1)]
            for m in range(1, p + 1):
                a.append(a[-1] * (-(p - m + 1)))   # a_m = a_{m-1} * -(p-m+1)
            for m in range(0, p + 1):
                pm = p - m
                q = [mp.mpc(0)] * (Nq + 1)
                am = a[m]
                hit = False
                for n in range(1, Nq + 1):
                    sn = s[n]
                    if sn != 0:
                        q[n] = sn * am / (mp.mpf(n) ** (m + 1))
                        hit = True
                if hit:
                    add(pm, q)
            # ---- n == 0 part: raise L by 1 ----
            if s[0] != 0:
                qr = [mp.mpc(0)] * (Nq + 1)
                qr[0] = s[0] / (p + 1)
                add(p + 1, qr)
        return LSeries(out, Nq)

    def value_L0(self, qC):
        r"""Evaluate the L^0 (finite, shuffle-regularised) part at qbar=qC."""
        s = self.terms.get(0, [mp.mpc(0)] * (self.Nq + 1))
        val = mp.mpc(0)
        qn = mp.mpc(1)
        for n in range(0, self.Nq + 1):
            val += s[n] * qn
            qn *= qC
        return val

    def value_at_L(self, qC, Lval):
        r"""Evaluate the FULL L-polynomial at qbar=qC and L=Lval, where L is the
        regularised value of the divergent F(omega_0)=int_{i inf}^tau 2 pi i dtau',
        i.e. L = log(qbar) = 2 pi i tau (tangential base point i inf -> drop log q_base).
        This is the AGW convention: divergent words evaluate to genuine log(qbar) powers,
        and the physical master keeps them (they combine into log-letters at the cusp)."""
        tot = mp.mpc(0)
        for p, s in self.terms.items():
            val = mp.mpc(0)
            qn = mp.mpc(1)
            for n in range(0, self.Nq + 1):
                val += s[n] * qn
                qn *= qC
            tot += val * (Lval ** p)
        return tot


# ---------------------------------------------------------------------------
# AGW form -> qbar-series C(qbar)  (the d qbar/qbar coefficient)
# ---------------------------------------------------------------------------
def series_for_form(name, fr):
    r"""C(qbar) power-series list for an AGW form named like 'omega_0_tau',
    'omega_2_z1_tau', 'omega_3_z3_2tau', 'eta_2_tau'.  C = (dtau-coeff)/(2 pi i)."""
    Nq = fr.Nq
    tpi = _twopii()
    if name == "omega_0_tau":
        C = [mp.mpc(0)] * (Nq + 1)
        C[0] = mp.mpc(1)        # omega_0 = 2 pi i d tau = dq/q  => C=1
        return C
    if name == "eta_2_tau":
        b2 = fr.b2()
        # eta_2 = b2 * dtau/(2 pi i)  => dtau-coeff = b2/(2 pi i)
        # C = dtau-coeff/(2 pi i) = b2/(2 pi i)^2
        return [x / (tpi * tpi) for x in b2]
    m = re.match(r"omega_(\d)_z(\d)_(2?)tau$", name)
    if not m:
        raise ValueError(f"unknown form {name!r}")
    k = int(m.group(1)); zi = int(m.group(2)); N = 2 if m.group(3) == "2" else 1
    if k == 2:
        g = fr.g2(zi, N)
        # dtau-coeff = N g^(2)(z,N tau)/(2 pi i)  => C = N g^(2)/(2 pi i)^2
        return [N * x / (tpi * tpi) for x in g]
    if k == 3:
        g = fr.g3(zi, N)
        # omega_3 dtau-coeff = (2pi)^{-1} * 2 g^(3)(z,N tau) * N/(2 pi i)
        # C = dtau-coeff/(2 pi i) = (2pi)^{-1} 2 N g^(3)/(2 pi i)^2
        pre = (2 * mp.pi) ** (-1) * 2 * N / (tpi * tpi)
        return [pre * x for x in g]
    raise ValueError(f"unsupported omega_{k}")


# ---------------------------------------------------------------------------
# Iterated integral of a single word and of a coefficient*word sum
# ---------------------------------------------------------------------------
def F_word(forms, fr):
    r"""F(omega_1,...,omega_k) = int omega_1 int omega_2 ... (eq.1494), via the
    recursion d/dq F(w1,w2,...) = C_{w1}(q)/q * F(w2,...).  Build inner-first:
    start S = 1, then for each form from the LAST to the FIRST:
       S <- regint( C_form * S ).
    Return the L^0 finite part at qbar = fr.qC."""
    Nq = fr.Nq
    S = LSeries.one(Nq)
    for name in reversed(forms):
        C = series_for_form(name, fr)
        S = S.mul_qseries(C).regint()
    return S.value_L0(fr.qC)


def F_word_L(forms, fr, Lval=None):
    r"""Like F_word but evaluating the FULL L-polynomial at L=log(qbar) (AGW convention,
    divergent omega_0/omega_2/eta_2 words kept as genuine log(qbar) powers).
    Lval defaults to log(qbar) = 2 pi i tau, here = mp.log(qC) on the principal branch;
    for the Euclidean negative-real qC use the explicit 2 pi i tau_C if provided."""
    Nq = fr.Nq
    S = LSeries.one(Nq)
    for name in reversed(forms):
        C = series_for_form(name, fr)
        S = S.mul_qseries(C).regint()
    if Lval is None:
        Lval = mp.log(fr.qC)
    return S.value_at_L(fr.qC, Lval)


# ---------------------------------------------------------------------------
# Parse an .mpl iter_int polynomial
# ---------------------------------------------------------------------------
def parse_mpl_poly(text):
    r"""Parse a Maple polynomial in iter_int([...]) words plus bare symbols
    (C_4_2, C_4_3, L1, L2) into a list of (coeff_str, token) where token is either
    a list of form-names (an iterated-integral word) or a ('SYM', name) /
    ('SYMxWORD', name, [forms]) for symbol-times-word products.

    Returns a list of dicts: {coeff: mpf-rational, kind, forms, sym}."""
    # normalise whitespace and stray backslashes
    t = text.strip()
    if t.startswith("(") and t.endswith(")"):
        t = t[1:-1]
    t = t.replace("\\", "").replace("\n", "").replace(" ", "")
    # split into additive terms at top-level +/-
    terms = []
    depth = 0
    cur = ""
    for ch in t:
        if ch == "[":
            depth += 1; cur += ch
        elif ch == "]":
            depth -= 1; cur += ch
        elif ch in "+-" and depth == 0 and cur != "" and cur[-1] not in "(*/^,":
            terms.append(cur); cur = ch
        else:
            cur += ch
    if cur:
        terms.append(cur)
    out = []
    for term in terms:
        term = term.strip()
        if not term:
            continue
        sign = 1
        if term[0] == "+":
            term = term[1:]
        elif term[0] == "-":
            sign = -1; term = term[1:]
        # extract coefficient (leading rational a/b or integer) and the symbolic body
        # patterns:  6*C_4_2*iter_int([...])  | C_4_3 | 1/2*C_4_2*iter_int([..]) | 24*iter_int([..])
        # split on '*'
        factors = term.split("*")
        coeff = mp.mpf(sign)
        sym = None
        forms = []
        i = 0
        while i < len(factors):
            f = factors[i]
            if f == "":
                i += 1; continue
            if re.match(r"^\d+/\d+$", f):
                num, den = f.split("/"); coeff *= mp.mpf(num) / mp.mpf(den)
            elif re.match(r"^\d+$", f):
                coeff *= mp.mpf(f)
            elif f in ("C_4_2", "C_4_3", "L1", "L2", "zeta_2", "zeta_3"):
                sym = f
            elif f.startswith("iter_int"):
                # f = iter_int([a,b,c]) possibly split across '*'? no, brackets keep it.
                inner = f[f.index("[") + 1:f.rindex("]")]
                forms = [w.strip() for w in inner.split(",") if w.strip()]
            else:
                raise ValueError(f"unparsed factor {f!r} in term {term!r}")
            i += 1
        out.append({"coeff": coeff, "sym": sym, "forms": forms})
    return out

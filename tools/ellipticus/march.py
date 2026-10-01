"""ellipticus.march — certified adaptive-Taylor transport of extended
inhomogeneous-PF systems with dlog-letter states (the workhorse,
engine i).

The transport core has been gated at 133d against independent quadrature;
this generalized form:
  * System takes the polyform dict (from ellipticus.extend.polyform) directly
    — no cache files, no frame names;
  * letters (a, i): appended states L' = Y_i/(z-a), a exact rational or an
    mp value produced from exact data at working precision;
  * transport: straight segments with |h| <= ratio * dist(exact pole set +
    letter poles), per-step tail assert < 1e-(wdps-8) (fail-closed), step cap;
    complex WAYPOINTS supported by chaining segments (the detour-tent
    discipline: route over a complex tent instead of through a
    near-singular stretch);
  * reg helpers: tangential-base-point (shuffle) regularization at an
    ANALYTIC base — J(z) = int [K-K(0)] dz'/z' + K(0) ln z — and the
    regular letter int K dz'/(z'-r) series base.

Measured gates: 133.3d at dps60 / 172.4d at dps100, two-precision stable.
"""
import mpmath as mp
from fractions import Fraction

from .curve import exact, _mpf_frac


def _coeffs_mp(strlist):
    return [mp.mpf(Fraction(c).numerator) / mp.mpf(Fraction(c).denominator)
            for c in strlist]


def _shift_desc(coeffs_desc, zc):
    """descending coeffs of p(z) -> ascending coeffs of p(zc+t) (Horner)."""
    out = [coeffs_desc[0]]
    for c in coeffs_desc[1:]:
        new = [None] * (len(out) + 1)
        new[0] = out[0] * zc + c
        for k in range(1, len(out)):
            new[k] = out[k] * zc + out[k - 1]
        new[len(out)] = out[-1]
        out = new
    return out


class System:
    """Numeric transportable system at current mp working precision.
    pf: dict from ellipticus.extend.polyform (n, Den, NumM, den_factors);
    letters: list of (a, i) — appended states L' = Y_i/(z-a)."""

    def __init__(self, pf, letters=(), dps=60):
        self.n = pf['n']
        self.Den = _coeffs_mp(pf['Den'])
        self.NumM = [[(_coeffs_mp(e) if e is not None else None)
                      for e in row] for row in pf['NumM']]
        self.letters = list(letters)
        self.poles = self._poles(pf, dps)

    @staticmethod
    def _poles(pf, dps):
        import sympy as sp
        out = []
        with mp.workdps(dps + 40):
            for f in pf['den_factors']:
                e = sp.sympify(f)
                syms = sorted(e.free_symbols, key=str)
                if not syms:
                    continue
                assert len(syms) == 1, f'multi-symbol den factor {f!r}'
                p = sp.Poly(e, syms[0])
                cs = [mp.mpf(sp.Rational(c).p) / mp.mpf(sp.Rational(c).q)
                      for c in p.all_coeffs()]
                if len(cs) > 1:
                    out += [mp.mpc(r) for r in
                            mp.polyroots(cs, maxsteps=400, extraprec=200)]
        return out

    def dist_to_sing(self, zc):
        ds = [abs(zc - p) for p in self.poles]
        ds += [abs(zc - a) for a, _ in self.letters]
        return min(ds) if ds else mp.mpf(1)

    def taylor_step(self, zc, Yfull, h, NT):
        """one Taylor step: returns (Y(zc+h), tailmag)."""
        n, nl = self.n, len(self.letters)
        Dsh = _shift_desc(self.Den, zc)
        degD = len(Dsh) - 1
        Ash = [[(_shift_desc(e, zc) if e is not None else None)
                for e in row] for row in self.NumM]
        degA = max(len(e) - 1 for row in Ash for e in row if e is not None)
        lets = [(zc - a) for a, _ in self.letters]
        Yc = [list(Yfull[:n])]
        Lc = [list(Yfull[n:])]
        D0 = Dsh[0]
        for m in range(NT):
            rhs = [mp.mpc(0)] * n
            for j in range(0, min(m, degA) + 1):
                Ym = Yc[m - j]
                for i in range(n):
                    row = Ash[i]
                    s = rhs[i]
                    for k in range(n):
                        e = row[k]
                        if e is not None and j < len(e):
                            s += e[j] * Ym[k]
                    rhs[i] = s
            for j in range(1, min(m + 1, degD) + 1):
                cf = Dsh[j] * (m + 1 - j)
                Ymj = Yc[m + 1 - j]
                for i in range(n):
                    rhs[i] -= cf * Ymj[i]
            inv = 1 / (D0 * (m + 1))
            Yc.append([rhs[i] * inv for i in range(n)])
            Ln = []
            for li, (a, isrc) in enumerate(self.letters):
                # source may be a Y component (isrc < n) or an earlier
                # letter state (n <= isrc < n + li): nested dlog words.
                if isrc < n:
                    src = Yc[m][isrc]
                else:
                    j = isrc - n
                    assert j < li, 'letter source must precede the letter'
                    src = Lc[m][j]
                Ln.append((src - m * Lc[m][li]) / (lets[li] * (m + 1)))
            Lc.append(Ln)
        allc = [Yc[m] + Lc[m] for m in range(NT + 1)]
        tot = list(allc[NT])
        for m in range(NT - 1, -1, -1):
            tot = [tot[i] * h + allc[m][i] for i in range(n + nl)]
        scale = max(max(abs(v) for v in tot), mp.mpf(10) ** (-30))
        tail = max(sum(abs(allc[NT - q][i]) * abs(h) ** (NT - q)
                       for q in range(3))
                   for i in range(n + nl)) / scale
        return tot, tail

    def transport(self, z_from, Yfull, z_to, NT=None, ratio=None,
                  tail_eps=None, maxstep=None):
        """straight-segment transport z_from -> z_to (either complex);
        complex waypoints (detours) = repeated calls.  Fail-closed tail
        assert per step; step-count fuse."""
        prec_dec = mp.mp.dps
        if NT is None:
            NT = int(prec_dec * 1.05) + 30
        if ratio is None:
            ratio = mp.mpf('0.10')
        if tail_eps is None:
            tail_eps = mp.mpf(10) ** (-(prec_dec - 8))
        zc = mp.mpc(z_from)
        zt = mp.mpc(z_to)
        Y = list(Yfull)
        nstep = 0
        while abs(zt - zc) > 0:
            dseg = zt - zc
            dist = self.dist_to_sing(zc)
            hmax = ratio * dist
            if maxstep is not None:
                hmax = min(hmax, mp.mpf(maxstep))
            h = dseg if abs(dseg) <= hmax else dseg / abs(dseg) * hmax
            Y, tail = self.taylor_step(zc, Y, h, NT)
            assert tail < tail_eps, \
                (f'transport tail {mp.nstr(tail, 5)} at z={mp.nstr(zc, 8)} '
                 f'(h={mp.nstr(h, 5)})')
            zc = zc + h
            nstep += 1
            assert nstep < 3000, 'step explosion'
        return Y, nstep

    def transport_path(self, z_from, Yfull, waypoints):
        """chained transport through a list of (possibly complex) waypoints
        — the detour-tent discipline."""
        zc = mp.mpc(z_from)
        Y = list(Yfull)
        total = 0
        for zt in waypoints:
            Y, ns = self.transport(zc, Y, zt)
            zc = mp.mpc(zt)
            total += ns
        return Y, total


# ------------------------------------------------- reg / letter series bases
def letter_base_log(K0, series, z1, nmax):
    """reg int_0^z1 K dz/z = K(0) ln z1 + sum_{n>=1} c_n z1^n / n
    (tangential-base-point / shuffle convention, unit tangent)."""
    v = K0 * mp.log(z1)
    for n in range(1, nmax + 1):
        v += series[n] * z1 ** n / n
    return v


def letter_base_r(series, z1, r, nmax):
    """int_0^z1 K dz/(z-r), r off (0, z1]:  sum_n c_n E_n,
    E_n = z1^n/n + r E_{n-1}, E_0 = ln((r-z1)/r)."""
    E = mp.log((r - z1) / r)
    v = series[0] * E
    for n in range(1, nmax + 1):
        E = z1 ** n / n + r * E
        v += series[n] * E
    return v


def seed_from_series(series_dict, z1, nmax):
    """kernel seeds K_k(z1) = sum_n c_n z1^n from exact local series."""
    return [sum(series_dict[k][n] * z1 ** n for n in range(nmax + 1))
            for k in sorted(series_dict)]

r"""membound.core — generic omega->0 memory-boundary omega-core engine.

Generalizes the single-constant c_M computation of the 5PM-2SF memory sector
(c_M = 1, reproduced here as the validation gate) to arbitrary 2- and
3-frequency memory omega-cores:

    J(eps) = (1/(2pi)^n) Int dw1..dwn  prod_ret (0^+ - i s_i w_i)^{p_i - q_i eps}
                                       prod_K   |w_j|^{-eps} K_eps(|w_j|)
                                       prod_poly w_i^{n_i}
  with the coupler frequency w_{n+1} = w_1 + ... + w_n (w3 = w1 + w2 for
  n_freq = 2; w4 = w1 + w2 + w3 for n_freq = 3), expanded analytically in eps
  to O(eps^2) BEFORE integration (avoids fixed-eps catastrophic cancellation),
  with |w|^{-eps} layers (the |w_i|^{+eps} in arXiv:2601.16256 eq. (27) is a
  sign typo; -eps is what reproduces c_M = 1).

Each frequency factor is split into its leading (eps^0) part, built with the
retarded sign structure (-i)^{p} w^{p} for EVERY flag (base_phase below), and
its order-eps logarithm
    log(0^+ - i s w) -> log|w| + i phi(sign w),
whose phase is set per factor by the flag presc in {'ret', 'adv', 'fey'}:
    phi_ret(s) = -(pi/2) s,   phi_adv(s) = +(pi/2) s,   phi_fey(s) = -(pi/2).
'ret' on every factor is the retarded assignment and is exact; 'adv' on every
factor is exact too, since w_i -> -w_i maps retarded to advanced and leaves
the integrals unchanged. 'fey' is a PHASE-ONLY SUBSTITUTION: the constant
Feynman order-eps phase on top of the retarded leading factor. It is a
diagnostic of how an assembled constant responds to the phase layer B, NOT
the Feynman propagator prescription, which changes the whole power to
(0^+ - i|w|)^{p - q eps} (leading factor -w1 w2 -> -|w1 w2| in I1^(M)); under
that assignment the leading cores change and the eps^-4 pole of I2^(M) does
not cancel, a case this engine does not implement. Mixed 'ret'/'adv' flags
likewise miss the overall sign (-1)^p of an advanced leading factor; use
uniform flags.
Per sign-sector the eps-expansion weight is derived AUTOMATICALLY:
    weight_k = proj[ base_phase * prod_i s_i^{n_i} * c_k ],
    c_0 = 1,  c_1 = A + iB,  c_2 = (A+iB)^2/2 + R,
    A = -sum_ret q_i log|w_i| - sum_K log|w_j|      (real log layer)
    B = -sum_ret q_i phi_i(s_i)                     (order-eps phase layer, per flag)
    R = sum_K kappa2(|w_j|)/K0(|w_j|)               (Bessel nu^2 layer)
    base_phase = (-i)^{sum_ret p_i}                 (retarded leading phase, used for every flag)
This reproduces the single-constant prototype's hard-coded
B1 = -(3 l1 + 3 l2 + l3) and the -2 pi^2 [s1==s2] term of B2 as special
cases (verified in the c_M gate).

proj in {'re','im'}: 're' = conservative (real-even) projection; 'im' is the
odd-in-v dissipative slot (NOT yet validated against anything — no published
values exist).

Cost law: quadrature = mp.quad tanh-sinh on panel breakpoints. kappa2 =
4th-order even FD with a dps-SCALED guard (no fixed precision guards; see
tools/baller, hygiene member dps_lint): workdps W = max(30, ceil(1.5*(dps+8))), h = 10^(-W/6)
=> kappa2 error ~ 10^(-2W/3) ~ 10^(-(dps+8)).
"""
import functools
import mpmath as mp

# Order-eps phase of log(0^+ - i s w) per flag; the ONLY place the flag enters.
# 'fey' swaps in the constant Feynman phase but NOT the Feynman leading factor:
# a phase-only diagnostic, not the Feynman propagator prescription (see the
# module docstring). The key name is kept for interface stability.
PHI = {
    'ret': lambda s: -mp.pi / 2 * s,
    'adv': lambda s: +mp.pi / 2 * s,
    'fey': lambda s: -mp.pi / 2,   # phase-only substitution, see above
}


class KCache:
    """K0 and kappa2 = (1/2) d^2/dnu^2 K_nu|_0 with dps-scaled FD guards."""

    def __init__(self, dps):
        self.dps = dps
        self.W = max(30, int(1.5 * (dps + 8)) + 1)
        self.h = mp.mpf(10) ** (-(self.W // 6))
        self.K0 = functools.lru_cache(maxsize=1_000_000)(self._k0)
        self.kap2 = functools.lru_cache(maxsize=1_000_000)(self._kap2)

    def _k0(self, x):
        return mp.besselk(0, x)

    def _kap2(self, x):
        # kap2 ~ (16 K_h - K_2h - 15 K_0)/(12 h^2), err O(h^4)  [prototype law]
        with mp.workdps(self.W):
            xx = mp.mpf(x)
            k0 = mp.besselk(0, xx)
            kh = mp.besselk(self.h, xx)
            k2h = mp.besselk(2 * self.h, xx)
            v = (16 * kh - k2h - 15 * k0) / (12 * self.h * self.h)
        return +v  # downcast to working dps


class OmegaCoreSpec:
    """Two- or three-frequency omega-core.

    n_freq: number of independent frequencies (2 or 3); the coupler index
            n_freq + 1 carries w_{n_freq+1} = w_1 + ... + w_{n_freq}
    ret   : list of (idx, p, q, presc) for (0^+ - i w_idx)^(p - q*eps),
            idx in 1..n_freq+1; presc in {'ret','adv','fey'} sets the
            order-eps phase only ('fey' = phase-only diagnostic, not a
            Feynman propagator; see module docstring)
    klines: list of idx carrying |w_idx|^{-eps} K_eps(|w_idx|)
    poly  : dict idx -> extra integer power of w_idx
    proj  : 're' | 'im'
    """

    FLAGS = ('ret', 'adv', 'fey')

    def __init__(self, name, ret, klines, poly=None, proj='re', n_freq=2):
        if n_freq not in (2, 3):
            raise ValueError(f"n_freq must be 2 or 3, got {n_freq}")
        if proj not in ('re', 'im'):
            raise ValueError(f"proj must be 're' or 'im', got {proj!r}")
        self.name = name
        # p is the integer power of the leading factor; q the eps-weight
        self.ret = [(int(i), int(p), q, str(f)) for i, p, q, f in ret]
        self.klines = list(klines)
        self.poly = {int(k): int(v) for k, v in dict(poly or {}).items()}
        self.proj = proj
        self.n_freq = n_freq
        self.idx = tuple(range(1, n_freq + 2))
        for idx, p, q, presc in self.ret:
            if isinstance(q, bool) or not isinstance(q, (int, float, mp.mpf)):
                raise ValueError(f"ret eps-weight q must be a number, got {q!r}")
            if presc not in self.FLAGS:
                raise ValueError(f"ret flag must be one of {self.FLAGS}, "
                                 f"got {presc!r}")
            if idx not in self.idx:
                raise ValueError(f"ret index {idx} outside 1..{n_freq + 1}")
        for idx in list(self.klines) + list(self.poly):
            if idx not in self.idx:
                raise ValueError(f"line index {idx} outside 1..{n_freq + 1}")
        # net integer power of each frequency (sign-carrying)
        self.npow = {i: 0 for i in self.idx}
        for idx, p, q, presc in self.ret:
            self.npow[idx] += p
        for idx, n in self.poly.items():
            self.npow[idx] += n
        self.base_phase = (-1j) ** sum(p for _, p, _, _ in self.ret)

    @property
    def layer_k(self):
        """|omega|^{-k eps} layer index = total frequency eps-weight."""
        return sum(q for _, _, q, _ in self.ret) + len(self.klines)

    # --- JSON round trip (the `python3 -m membound --spec core.json` schema) --
    def to_dict(self):
        """Plain-JSON form: {name, n_freq, ret:[[idx,p,q,flag],...], klines,
        poly:{str(idx): n}, proj}."""
        return {"name": self.name, "n_freq": self.n_freq,
                "ret": [list(r) for r in self.ret],
                "klines": list(self.klines),
                "poly": {str(k): v for k, v in self.poly.items()},
                "proj": self.proj}

    @classmethod
    def from_dict(cls, d):
        """Inverse of to_dict. `poly` may be a {idx: n} mapping (string or
        integer keys), a list of [idx, n] pairs, or absent/null."""
        unknown = set(d) - {"name", "n_freq", "ret", "klines", "poly", "proj",
                            "comment", "_comment"}
        if unknown:
            raise ValueError(f"unknown core-spec keys: {sorted(unknown)}")
        for key in ("name", "ret", "klines"):
            if key not in d:
                raise ValueError(f"core spec is missing required key {key!r}")
        poly = d.get("poly") or {}
        if isinstance(poly, list):
            poly = {int(i): int(n) for i, n in poly}
        return cls(name=str(d["name"]),
                   ret=[tuple(r) for r in d["ret"]],
                   klines=[int(j) for j in d["klines"]],
                   poly=poly,
                   proj=d.get("proj", "re"),
                   n_freq=int(d.get("n_freq", 2)))


# --- sign-sector domains, n_freq = 2 -----------------------------------------
# Domain = (tag, s1, s2). Coordinates (x, y) in (0, inf)^2 map to
# (|w1|, |w2|, |w3|, s3):
#   same  (s1==s2): |w1|=x, |w2|=y, |w3|=x+y, s3=s1
#   oppgt (s1!=s2, |w1|>|w2|): |w2|=x, |w3|=y, |w1|=x+y, s3=s1
#   opplt (s1!=s2, |w1|<|w2|): |w1|=x, |w3|=y, |w2|=x+y, s3=s2
def domains():
    out = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            if s1 == s2:
                out.append(('same', s1, s2))
            else:
                out.append(('oppgt', s1, s2))
                out.append(('opplt', s1, s2))
    return out


def _coords(tag, s1, s2, x, y):
    if tag == 'same':
        return x, y, x + y, s1
    if tag == 'oppgt':
        return x + y, x, y, s1
    return x, x + y, y, s2


# --- sign-sector domains, n_freq = 3 (sign octants) --------------------------
# Domain = (tag, s1, s2, s3) over the 8 sign octants of (w1, w2, w3), coupler
# w4 = w1 + w2 + w3. Writing a_i = |w_i|:
#   same (s1==s2==s3): a1=x, a2=y, a3=z; a4 = x+y+z, s4 = s1.
#   A mixed octant has a lone index l (the minority sign) and a pair {i, j},
#   i < j; it splits at the |w4| = 0 locus:
#     gt (a_i + a_j > a_l): a_l = x, a4 = y, a_i = (x+y) u, a_j = (x+y)(1-u)
#        with u = z in (0, 1); s4 = s_pair; Jacobian (x + y).
#     lt (a_i + a_j < a_l): a_i = x, a_j = y, a4 = z, a_l = x+y+z; s4 = s_l;
#        Jacobian 1.
# 2 uniform + 6 mixed x 2 = 14 domains. On both mixed pieces the coupler
# magnitude a4 is an explicit coordinate, so its singular edge (K line or
# negative power) sits on a panel boundary where tanh-sinh clusters nodes.
def domains3():
    out = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            for s3 in (1, -1):
                if s1 == s2 == s3:
                    out.append(('same', s1, s2, s3))
                else:
                    out.append(('gt', s1, s2, s3))
                    out.append(('lt', s1, s2, s3))
    return out


def _lone_pair(s1, s2, s3):
    """Mixed-octant split: 1-based (lone, pair_i, pair_j), pair_i < pair_j."""
    s = (s1, s2, s3)
    for k in (0, 1, 2):
        i, j = [m for m in (0, 1, 2) if m != k]
        if s[i] == s[j] != s[k]:
            return k + 1, i + 1, j + 1
    raise ValueError(f"not a mixed octant: {(s1, s2, s3)}")


def _coords3(tag, s1, s2, s3, x, y, z):
    """-> (a1, a2, a3, a4, s4, jacobian)."""
    if tag == 'same':
        return x, y, z, x + y + z, s1, mp.mpf(1)
    lone, i, j = _lone_pair(s1, s2, s3)
    s = {1: s1, 2: s2, 3: s3}
    a = {}
    if tag == 'gt':
        a[lone] = x
        a[i] = (x + y) * z
        a[j] = (x + y) * (1 - z)
        return a[1], a[2], a[3], y, s[i], x + y
    # 'lt'
    a[i] = x
    a[j] = y
    a[lone] = x + y + z
    return a[1], a[2], a[3], z, s[lone], mp.mpf(1)


def _panel_sig(dom):
    """Quadrature-panel signature: the mixed 'gt' pieces of a 3-frequency
    core integrate their third axis over (0, 1); everything else uses the
    magnitude panels on every axis."""
    return 'simplex' if len(dom) == 4 and dom[0] == 'gt' else 'prod'


def _weight(spec, absw, sgn, base, order, kc, proj):
    """Apply the automatic eps-layer weight c_k to a sector base value."""
    base_phase = spec.base_phase
    if order == 0:
        return proj(base_phase) * base
    A = mp.mpf(0)
    B = mp.mpf(0)
    for idx, p, q, presc in spec.ret:
        A -= q * mp.log(absw[idx])
        B -= q * PHI[presc](sgn[idx])
    for j in spec.klines:
        A -= mp.log(absw[j])
    if order == 1:
        return proj(base_phase * mp.mpc(A, B)) * base
    # order 2: (A+iB)^2/2 + R
    R = mp.mpf(0)
    for j in spec.klines:
        R += kc.kap2(absw[j]) / kc.K0(absw[j])
    c2 = mp.mpc(A, B) ** 2 / 2 + R
    return proj(base_phase * c2) * base


def _sector_base(spec, absw, sgn, kc):
    """|w|^n powers, sector signs, K0 lines (the order-0 sector integrand)."""
    mag = mp.mpf(1)
    ssign = 1
    for i in spec.idx:
        n = spec.npow[i]
        if n:
            mag *= absw[i] ** n
            if n % 2:
                ssign *= sgn[i]
    kprod = mp.mpf(1)
    for j in spec.klines:
        kprod *= kc.K0(absw[j])
    return ssign * mag * kprod


def make_integrand(spec, dom, order, kc):
    """Scalar integrand for one domain and one eps-order (0, 1, 2)."""
    proj = (lambda z: mp.re(z)) if spec.proj == 're' else (lambda z: mp.im(z))

    if spec.n_freq == 3:
        tag, s1, s2, s3 = dom

        def f3(x, y, z):
            a1, a2, a3, a4, s4, jac = _coords3(tag, s1, s2, s3, x, y, z)
            absw = {1: a1, 2: a2, 3: a3, 4: a4}
            sgn = {1: s1, 2: s2, 3: s3, 4: s4}
            base = jac * _sector_base(spec, absw, sgn, kc)
            return _weight(spec, absw, sgn, base, order, kc, proj)

        return f3

    tag, s1, s2 = dom

    def f(x, y):
        a1, a2, a3, s3 = _coords(tag, s1, s2, x, y)
        absw = {1: a1, 2: a2, 3: a3}
        sgn = {1: s1, 2: s2, 3: s3}
        base = _sector_base(spec, absw, sgn, kc)
        return _weight(spec, absw, sgn, base, order, kc, proj)

    return f


_PROBES = [("0.7", "1.3"), ("0.31", "2.1"), ("1.9", "0.45")]
# 3-frequency probes: the third coordinate stays in (0, 1) so one probe set is
# valid both as a magnitude and as the simplex coordinate u of 'gt' pieces.
_PROBES3 = [("0.7", "1.3", "0.45"), ("0.31", "2.1", "0.6"),
            ("1.9", "0.45", "0.85")]


def domain_classes(spec, order, kc):
    """Group the sign-domains into equal-integrand classes by probing.

    Correctness-relevant grouping: two domains merge only if the FULL scalar
    integrand agrees at 3 probe points to (dps-3) digits AND they share a
    quadrature-panel signature (a 'gt' piece never merges with a
    product-coordinate domain — their third-axis panels differ). Returns
    [(representative_domain, multiplicity, member_list)].
    """
    if spec.n_freq == 3:
        doms, probes = domains3(), _PROBES3
    else:
        doms, probes = domains(), _PROBES
    keys = {}
    tol_digits = mp.mp.dps - 3
    for d in doms:
        f = make_integrand(spec, d, order, kc)
        key = (_panel_sig(d),) + tuple(
            mp.nstr(f(*[mp.mpf(p) for p in pt]), tol_digits) for pt in probes)
        keys.setdefault(key, []).append(d)
    return [(members[0], len(members), members) for members in keys.values()]


def compute_core(spec, dps, brk, maxdegree, log=print, kc=None,
                 orders=(0, 1, 2)):
    """J^(k) for k in orders (default 0, 1, 2): the eps-layer moments of the
    omega-core, i.e. the Laurent coefficients of J(eps) = sum_k J^(k) eps^k.

    Returns dict order -> mpf/mpc value INCLUDING the 1/(2pi)^n_freq measure.
    brk applies to every magnitude axis; the (0, 1) simplex axis of the
    3-frequency 'gt' pieces uses a single (0, 1) panel.
    """
    import time
    kc = kc or KCache(dps)
    inv = 1 / (2 * mp.pi) ** spec.n_freq
    ubrk = [mp.mpf(0), mp.mpf(1)]
    out = {}
    for order in orders:
        classes = domain_classes(spec, order, kc)
        tot = mp.mpf(0)
        for rep, mult, members in classes:
            t0 = time.time()
            f = make_integrand(spec, rep, order, kc)
            if spec.n_freq == 3:
                third = ubrk if rep[0] == 'gt' else brk
                val = mp.quad(f, brk, brk, third, maxdegree=maxdegree)
            else:
                val = mp.quad(f, brk, brk, maxdegree=maxdegree)
            dt = time.time() - t0
            # drop numerically-dead classes cleanly (im-part of real sector etc.)
            tot += mult * val
            log(f"    [{spec.name} k={order}] class rep={rep} x{mult} "
                f"quad={mp.nstr(val, min(20, dps))}  ({dt:.1f}s)")
        out[order] = inv * tot
    return out

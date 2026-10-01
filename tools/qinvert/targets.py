"""qinvert.targets — statistic definitions as first-class objects.

Each statistic exposes
  value(xs_sorted)  : exact evaluation on a sorted list of Fractions, and
  support(n)        : the order-statistic sensitivity structure — the 1-based
                      indices and rational weights the statistic depends on for
                      a dataset of size n (None for dense statistics like Mean).
value() is DERIVED from support() for every quantile, so the solver's
sensitivity structure and the evaluator can never disagree; the tests
cross-check both against independently written textbook formulas and
published R reference values.

Every accept/reject quantity is a fractions.Fraction.  No float participates
in any statistic value anywhere in this module (the publisher conventions
handled so far — the documented QNTLDEF quantile definitions applied to
decimal-string scores — have exact-rational semantics, so there are zero
documented float sites in qinvert).

Provenance
----------
- SAS QNTLDEF=1..5: Base SAS 9.4 Procedures Guide: Statistical Procedures,
  chapter "The UNIVARIATE Procedure", section "Calculating Percentiles"
  (QNTLDEF= / PCTLDEF= definitions 1-5; PROC MEANS shares the definitions,
  default QNTLDEF=5).  QNTLDEF=5 is additionally validated end-to-end
  against the CMS Tech Notes K-5/K-6 fence tables (Tech Notes 2024 final
  p.149-150).
- R types 1-9: Hyndman, R.J. & Fan, Y. (1996), "Sample Quantiles in
  Statistical Packages", The American Statistician 50(4), 361-365 (the R
  taxonomy; R's ?quantile documents the same nine types).  Documented
  equivalences: SAS1=R4, SAS2=R3, SAS3=R1, SAS4=R6, SAS5=R2 (H&F sec. 2 +
  SAS doc notes); asserted by the test battery.
- Tukey fence algebra: lo = Q1 - m*IQR = (1+m)Q1 - m*Q3,
  hi = Q3 + m*IQR = (1+m)Q3 - m*Q1; inversion
  Q1 = ((1+m)lo + m hi)/(1+2m), Q3 = (m lo + (1+m) hi)/(1+2m).
  m=3 specialization (4lo+3hi)/7, (3lo+4hi)/7 is validated by the
  planted-truth battery.
"""

from fractions import Fraction as F

__all__ = [
    "Interval", "Statistic", "Quantile", "OrderStat", "Mean", "TukeyFence",
    "LinearCombo", "RangeStat", "MinStat", "MaxStat",
    "invert_fence_pair", "as_fraction", "evaluate_targets", "targets_exact",
    "SAS_DEFINITIONS", "R_DEFINITIONS", "ALL_DEFINITIONS",
]

SAS_DEFINITIONS = ("sas1", "sas2", "sas3", "sas4", "sas5")
R_DEFINITIONS = tuple(f"r{i}" for i in range(1, 10))
ALL_DEFINITIONS = SAS_DEFINITIONS + R_DEFINITIONS


def as_fraction(x):
    """Exact conversion: Fraction/int pass through, decimal strings parsed
    exactly ('1.41' -> 141/100).  Floats are rejected — a float here would
    mean a convention we have not certified (no replayed pipeline is
    float-semantic; see module docstring)."""
    if isinstance(x, F):
        return x
    if isinstance(x, int):
        return F(x)
    if isinstance(x, str):
        return F(x)
    raise TypeError(f"refusing inexact type {type(x).__name__}: {x!r}")


class Interval:
    """Published interval target [lo, hi] (None = unbounded side).
    Used for capped / inequality published statistics."""

    def __init__(self, lo=None, hi=None):
        self.lo = as_fraction(lo) if lo is not None else None
        self.hi = as_fraction(hi) if hi is not None else None

    def contains(self, v):
        return ((self.lo is None or v >= self.lo) and
                (self.hi is None or v <= self.hi))

    def __repr__(self):
        return f"Interval({self.lo}, {self.hi})"


def _floor(x):
    """Exact floor of a Fraction."""
    return x.numerator // x.denominator


def _cl(i, n):
    """Clamp a 1-based index into [1, n]."""
    return min(max(i, 1), n)


def _support_sas(qdef, n, p):
    """Order-stat support of SAS QNTLDEF=qdef at size n, probability p.

    Notation (SAS 9.4 Proc Guide, 'Calculating Percentiles'): np = n*p,
    j = floor(np), g = np - j.  Returns ((idx, weight), ...) 1-based.
    """
    np_ = n * p
    j = _floor(np_)
    g = np_ - j
    if qdef == 1:
        # QNTLDEF=1: weighted average at x_{np}: y = (1-g)x_j + g x_{j+1},
        # with x_0 := x_1 and x_{n+1} := x_n (SAS doc edge convention).
        if g == 0:
            return ((_cl(j, n), F(1)),)
        jj, jj1 = _cl(j, n), _cl(j + 1, n)
        if jj == jj1:
            return ((jj, F(1)),)
        return ((jj, 1 - g), (jj1, g))
    if qdef == 2:
        # QNTLDEF=2: observation numbered closest to np; tie g=1/2 resolved
        # to x_j when j is even, else x_{j+1} (SAS doc).
        if g == F(1, 2):
            i = j if j % 2 == 0 else j + 1
        elif g < F(1, 2):
            i = j
        else:
            i = j + 1
        return ((_cl(i, n), F(1)),)
    if qdef == 3:
        # QNTLDEF=3: empirical distribution function: y = x_j if g == 0
        # else x_{j+1}.
        i = j if g == 0 else j + 1
        return ((_cl(i, n), F(1)),)
    if qdef == 4:
        # QNTLDEF=4: weighted average aimed at x_{(n+1)p}: h = (n+1)p,
        # j = floor(h), g = h - j; y = (1-g)x_j + g x_{j+1}, clamped.
        h = (n + 1) * p
        j4 = _floor(h)
        g4 = h - j4
        if g4 == 0:
            return ((_cl(j4, n), F(1)),)
        jj, jj1 = _cl(j4, n), _cl(j4 + 1, n)
        if jj == jj1:
            return ((jj, F(1)),)
        return ((jj, 1 - g4), (jj1, g4))
    if qdef == 5:
        # QNTLDEF=5 (SAS default): EDF with averaging: y = (x_j + x_{j+1})/2
        # if g == 0 else x_{j+1}.  (module docstring provenance).
        if g == 0:
            jj, jj1 = _cl(j, n), _cl(j + 1, n)
            if jj == jj1:
                return ((jj, F(1)),)
            return ((jj, F(1, 2)), (jj1, F(1, 2)))
        return ((_cl(j + 1, n), F(1)),)
    raise ValueError(f"unknown SAS QNTLDEF {qdef}")


def _support_r(rtype, n, p):
    """Order-stat support of Hyndman-Fan type rtype (1-9) at size n, prob p.

    Discontinuous types 1-3: j = floor(np + m), g = np + m - j with
    m = 0 (types 1,2) or -1/2 (type 3); Q per-type step rule.
    Continuous types 4-9: h per type, j = floor(h), g = h - j,
    Q = (1-g)x_j + g x_{j+1}, clamped to [x_1, x_n] outside [1, n].
    (Hyndman & Fan 1996, eqs. (1)-(9) / R ?quantile.)
    """
    if rtype in (1, 2, 3):
        m = F(0) if rtype in (1, 2) else F(-1, 2)
        h = n * p + m
        j = _floor(h)
        g = h - j
        if rtype == 1:
            i = j if g == 0 else j + 1
            return ((_cl(i, n), F(1)),)
        if rtype == 2:
            if g == 0:
                jj, jj1 = _cl(j, n), _cl(j + 1, n)
                if jj == jj1:
                    return ((jj, F(1)),)
                return ((jj, F(1, 2)), (jj1, F(1, 2)))
            return ((_cl(j + 1, n), F(1)),)
        # type 3: nearest order statistic, ties to even j
        i = j if (g == 0 and j % 2 == 0) else j + 1
        return ((_cl(i, n), F(1)),)
    if rtype == 4:
        h = n * p
    elif rtype == 5:
        h = n * p + F(1, 2)
    elif rtype == 6:
        h = (n + 1) * p
    elif rtype == 7:
        h = (n - 1) * p + 1
    elif rtype == 8:
        h = (n + F(1, 3)) * p + F(1, 3)
    elif rtype == 9:
        h = (n + F(1, 4)) * p + F(3, 8)
    else:
        raise ValueError(f"unknown R type {rtype}")
    if h <= 1:
        return ((1, F(1)),)
    if h >= n:
        return ((n, F(1)),)
    j = _floor(h)
    g = h - j
    if g == 0:
        return ((j, F(1)),)
    return ((j, 1 - g), (j + 1, g))


class Statistic:
    """Base statistic.  dense=True means the value depends on every data
    point (e.g. Mean), so the order-stat placement engine cannot pin it and
    the solver treats it as a sum constraint / verification target."""

    dense = False

    def value(self, xs_sorted):
        raise NotImplementedError

    def support(self, n):
        """((idx, weight), ...) 1-based, or None for dense statistics."""
        return None


class Quantile(Statistic):
    """Quantile of probability p under a named definition.

    definition in SAS_DEFINITIONS ('sas1'..'sas5', SAS QNTLDEF=1..5) or
    R_DEFINITIONS ('r1'..'r9', Hyndman-Fan types 1-9).  See module docstring
    for provenance.
    """

    def __init__(self, p, definition="sas5"):
        self.p = as_fraction(p)
        if not (0 < self.p < 1):
            raise ValueError("p must be in (0,1)")
        definition = str(definition).lower()
        if definition not in ALL_DEFINITIONS:
            raise ValueError(f"unknown definition {definition!r}; "
                             f"choose from {ALL_DEFINITIONS}")
        self.definition = definition

    def support(self, n):
        if n < 1:
            raise ValueError("empty data")
        if self.definition.startswith("sas"):
            return _support_sas(int(self.definition[3:]), n, self.p)
        return _support_r(int(self.definition[1:]), n, self.p)

    def value(self, xs_sorted):
        sup = self.support(len(xs_sorted))
        return sum(w * xs_sorted[i - 1] for i, w in sup)

    def __repr__(self):
        return f"Quantile({self.p}, {self.definition!r})"


class OrderStat(Statistic):
    """The j-th order statistic (1-based).  Negative j counts from the top:
    j=-1 is the maximum.  MinStat()/MaxStat() are conveniences."""

    def __init__(self, j):
        if j == 0:
            raise ValueError("j is 1-based; 0 is invalid")
        self.j = j

    def index(self, n):
        i = self.j if self.j > 0 else n + 1 + self.j
        if not (1 <= i <= n):
            raise ValueError(f"index {self.j} out of range for n={n}")
        return i

    def support(self, n):
        return ((self.index(n), F(1)),)

    def value(self, xs_sorted):
        return xs_sorted[self.index(len(xs_sorted)) - 1]

    def __repr__(self):
        return f"OrderStat({self.j})"


def MinStat():
    return OrderStat(1)


def MaxStat():
    return OrderStat(-1)


class Mean(Statistic):
    """Arithmetic mean — dense: depends on every value, so the solver turns
    an exact Mean target into an exact sum constraint on the delta."""

    dense = True

    def value(self, xs_sorted):
        return F(sum(xs_sorted), len(xs_sorted))

    def __repr__(self):
        return "Mean()"


class TukeyFence(Statistic):
    """Tukey fence: side='lo' -> Q(p_lo) - mult*IQR, side='hi' ->
    Q(p_hi) + mult*IQR, where IQR = Q(p_hi) - Q(p_lo), quantiles under
    `definition`.  cap (optional) truncates the published value:
    lo-fence value = max(raw, cap); hi-fence value = min(raw, cap)
    (the CMS convention: percent measures capped into [0, 100], rates
    capped below at 0 — Tech Notes 2024 p.149-150).

    mult=3 is the OUTER fence; mult=3/2 the inner.
    """

    def __init__(self, side, mult=3, definition="sas5",
                 p_lo=F(1, 4), p_hi=F(3, 4), cap=None):
        if side not in ("lo", "hi"):
            raise ValueError("side must be 'lo' or 'hi'")
        self.side = side
        self.mult = as_fraction(mult)
        self.definition = definition
        self.p_lo = as_fraction(p_lo)
        self.p_hi = as_fraction(p_hi)
        self.q_lo = Quantile(self.p_lo, definition)
        self.q_hi = Quantile(self.p_hi, definition)
        self.cap = as_fraction(cap) if cap is not None else None

    def pair_key(self):
        """Fences with the same pair_key form an invertible lo/hi pair."""
        return (self.mult, self.definition, self.p_lo, self.p_hi)

    def raw_value(self, xs_sorted):
        q1 = self.q_lo.value(xs_sorted)
        q3 = self.q_hi.value(xs_sorted)
        m = self.mult
        # lo = Q1 - m(Q3-Q1) = (1+m)Q1 - mQ3 ; hi = Q3 + m(Q3-Q1) = (1+m)Q3 - mQ1
        if self.side == "lo":
            return (1 + m) * q1 - m * q3
        return (1 + m) * q3 - m * q1

    def value(self, xs_sorted):
        raw = self.raw_value(xs_sorted)
        if self.cap is None:
            return raw
        return max(raw, self.cap) if self.side == "lo" else min(raw, self.cap)

    def __repr__(self):
        return (f"TukeyFence({self.side!r}, mult={self.mult}, "
                f"definition={self.definition!r}, cap={self.cap})")


class LinearCombo(Statistic):
    """Linear combination sum(coeff * stat).  Dense if any term is dense.
    The solver can INVERT (drive enumeration with) a LinearCombo only when
    it has exactly two single-index order-stat terms (e.g. Range = max - min,
    2*midrange); other combos are verification-only targets (loud cap note).
    """

    def __init__(self, terms):
        self.terms = [(as_fraction(c), s) for c, s in terms]
        self.dense = any(s.dense for _, s in self.terms)

    def value(self, xs_sorted):
        return sum(c * s.value(xs_sorted) for c, s in self.terms)

    def __repr__(self):
        return f"LinearCombo({self.terms!r})"


def RangeStat():
    """Range = max - min (a 2-term order-stat LinearCombo the solver can
    invert by gridding the min)."""
    return LinearCombo([(F(1), MaxStat()), (F(-1), MinStat())])


def invert_fence_pair(pub_lo, pub_hi, mult, cap_lo, cap_hi):
    """Given published lo/hi fence values (Fractions or None), fence
    multiplier m, and the caps, return (q_lo_req, q_hi_req, mode):

      mode 'both'    : both fences uncapped -> both quantiles pinned:
                       Q1 = ((1+m)lo + m hi)/(1+2m), Q3 = (m lo + (1+m) hi)/(1+2m)
      mode 'hi_only' : only hi = (1+m)Q3 - m Q1 constrains (lo absent/capped)
      mode 'lo_only' : only lo = (1+m)Q1 - m Q3 constrains
      mode 'none'    : both sides capped/absent (inequalities only)

    A side is 'capped' when its published value EQUALS its cap: the raw
    fence then only satisfies an inequality, not an equation (validated by
    planted-truth batteries).
    """
    m = as_fraction(mult)
    lo_capped = (pub_lo is not None and cap_lo is not None and
                 as_fraction(pub_lo) == as_fraction(cap_lo))
    hi_capped = (pub_hi is not None and cap_hi is not None and
                 as_fraction(pub_hi) == as_fraction(cap_hi))
    if pub_lo is None or pub_hi is None:
        if pub_hi is not None and not hi_capped:
            return None, None, "hi_only"
        if pub_lo is not None and not lo_capped:
            return None, None, "lo_only"
        return None, None, "none"
    if not lo_capped and not hi_capped:
        lo, hi = as_fraction(pub_lo), as_fraction(pub_hi)
        d = 1 + 2 * m
        return ((1 + m) * lo + m * hi) / d, (m * lo + (1 + m) * hi) / d, "both"
    if lo_capped and not hi_capped:
        return None, None, "hi_only"
    if hi_capped and not lo_capped:
        return None, None, "lo_only"
    return None, None, "none"


def evaluate_targets(xs_sorted, targets):
    """[(value, ok), ...] parallel to targets; ok is the exact accept bit.
    targets: [(Statistic, published)] with published a Fraction-able exact
    value or an Interval.

    A statistic UNDEFINED on this array (e.g. an OrderStat whose index lies
    outside [1, n], or any statistic on empty data) evaluates to
    (None, False): it cannot match any published value, and evaluation must
    never crash a solve/stack run over one unsatisfiable dataset size
    (receipt: tests/test_all.py::TestOrderStatOutOfRange)."""
    out = []
    for stat, pub in targets:
        try:
            v = stat.value(xs_sorted)
        except ValueError:
            out.append((None, False))
            continue
        if isinstance(pub, Interval):
            out.append((v, pub.contains(v)))
        else:
            out.append((v, v == as_fraction(pub)))
    return out


def targets_exact(xs_sorted, targets):
    """True iff every target matches exactly.  This is the authoritative
    verification gate the solver applies to every witness."""
    return all(ok for _, ok in evaluate_targets(xs_sorted, targets))

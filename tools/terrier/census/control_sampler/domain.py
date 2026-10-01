"""Fundamental-domain abstraction for the control sampler.

The H0-a control is sampled by REJECTION on the SAME fundamental domain
as the census dedup. CONTRACT: the dedup module (Engines A+B) exports its
fundamental domain as an object with the interface below; the sampler
consumes THAT object unchanged. The sampler never constructs its own
domain for a census run -- the synthetic domains here exist only for the
fixed-seed test battery.

Interface (duck-typed):
  box()        -> (lo, hi): lists of Fractions, the bounding box
  contains(x)  -> bool, EXACT rational membership for rational x
"""
from fractions import Fraction


class BoxDomain:
    """A rectangular fundamental domain: membership == bounding box.
    Half-open [lo, hi) per axis (measure-zero convention only)."""

    def __init__(self, lo, hi):
        self.lo = [Fraction(v) for v in lo]
        self.hi = [Fraction(v) for v in hi]
        assert len(self.lo) == len(self.hi)
        assert all(a < b for a, b in zip(self.lo, self.hi))

    @property
    def dim(self):
        return len(self.lo)

    def box(self):
        return self.lo, self.hi

    def contains(self, x):
        return all(a <= v < b for v, a, b in zip(x, self.lo, self.hi))

    def volume(self):
        v = Fraction(1)
        for a, b in zip(self.lo, self.hi):
            v *= b - a
        return v


class PredicateDomain(BoxDomain):
    """Bounding box + exact-rational membership predicate (e.g. a modular
    fundamental domain handed over by the dedup engine)."""

    def __init__(self, lo, hi, predicate):
        super().__init__(lo, hi)
        self._pred = predicate

    def contains(self, x):
        return super().contains(x) and bool(self._pred(x))

"""baller.tripwire — the dual-path HALT-NOT-AVERAGE discipline as one
reusable component (the pattern behind cport's dual-run oracle, the
coalescer / nestor.dispersion two-method gates and wayfinder's digit gates).

Law: two independent routes compute the same quantity; if they agree to the
bar, the PRIMARY result is returned untouched; if they do not, computation
HALTS with DualPathDisagreement carrying both values — results are never
averaged, never "reconciled", never silently picked.
"""
import math
import zlib


class DualPathDisagreement(RuntimeError):
    """Named halt: the two routes disagreed below the bar (or a route went
    numerically meaningless — NaN/inf). Carries .primary, .oracle, .digits —
    do not use either value."""

    def __init__(self, msg, primary=None, oracle=None, digits=None):
        super().__init__(msg)
        self.primary, self.oracle, self.digits = primary, oracle, digits


def agree_digits(a, b, working_digits=None):
    """Significant digits of agreement between two REAL values, computed in
    mpmath at `working_digits` precision (default: enough to resolve the
    caller's bar — pass bar+20 from dual(); bare calls get 60). Guards:
    float64 collapse (two mpf disagreeing at digit 20 must never read as
    exact equality), NaN passing every bar, and tiny-vs-zero asymmetry (denominator
    is now max(|a|,|b|) — symmetric, so route order cannot decide whether the
    wire trips). Complex/ball inputs are REJECTED typed: compare mids/parts
    explicitly (documented boundary).
    Returns 9999.0 on exact equality; 0.0 when exactly one value is 0 (total
    relative disagreement — any bar >= 1 halts); -inf on non-finite values."""
    import mpmath
    # flint arb carries _mpf_, so mpmathify(arb) silently
    # yields the MIDPOINT — the ball's radius would be DISCARDED. Balls are
    # rejected typed: the caller compares mids explicitly (and owns the radius story).
    for v in (a, b):
        if hasattr(v, "rad") or hasattr(v, "mid"):
            raise ValueError(
                f"agree_digits: ball input ({type(v).__name__}) — a ball's "
                f"radius cannot ride through a scalar comparison; compare "
                f"mids explicitly and watch radii with RadiusWatch")
    wd = int(working_digits or 60)
    with mpmath.workdps(wd):
        try:
            ma, mb = mpmath.mpmathify(a), mpmath.mpmathify(b)
        except (TypeError, ValueError) as e:
            raise ValueError(
                f"agree_digits: inputs must be real scalars (got {type(a).__name__}"
                f"/{type(b).__name__}); compare ball mids / complex parts "
                f"explicitly") from e
        if isinstance(ma, mpmath.mpc) or isinstance(mb, mpmath.mpc):
            raise ValueError(
                "agree_digits: complex input — compare real/imag parts (or "
                "ball mids) explicitly; the tripwire bar is per-component")
        for v in (ma, mb):
            if not mpmath.isfinite(v):
                return -math.inf          # NaN/inf NEVER passes a bar
        if ma == mb:
            return 9999.0
        den = max(abs(ma), abs(mb))
        if den == 0:
            return 9999.0                 # both exactly zero
        return float(-mpmath.log10(abs(ma - mb) / den))


def dual(primary, oracle, digits, label="dual"):
    """Run both callables; halt unless they agree to `digits` significant
    digits (compared at digits+20 working precision — never through float64).
    Returns (primary_value, certificate_dict)."""
    p = primary()
    o = oracle()
    d = agree_digits(p, o, working_digits=digits + 20)
    if not (d >= digits):                 # NaN-safe: non-finite halts
        raise DualPathDisagreement(
            f"{label}: routes agree to only {d:.1f} digits (bar {digits}) — "
            f"HALT, never average (primary={p!r}, oracle={o!r})",
            primary=p, oracle=o, digits=d)
    return p, {"label": label, "agree_digits": min(d, 9999.0), "bar": digits}


class Tripwire:
    """Sampled dual-running for hot loops (the cport 1-in-N pattern): a
    deterministic crc32 selector dual-runs the oracle on a fixed fraction of
    calls; any disagreement below the bar halts the whole computation."""

    def __init__(self, oracle, digits, every=15, label="tripwire"):
        self.oracle, self.digits, self.every, self.label = oracle, digits, every, label
        self.calls = 0
        self.checked = 0
        self.worst = 9999.0

    def selected(self, key):
        """Deterministic selection: crc32 of the key's repr, 1-in-`every`.
        An unreprable arg falls back to (type-name, call-count) — the
        wire must never crash on the caller's objects."""
        try:
            blob = repr(key)
        except Exception:
            blob = f"unreprable:{type(key).__name__}:{self.calls}"
        return zlib.crc32(blob.encode()) % self.every == 0

    def __call__(self, primary_fn, *args, **kw):
        self.calls += 1
        p = primary_fn(*args, **kw)
        # the call counter is ALWAYS in the key — any constant
        # call signature (fixed config arg, closure-carried state) would
        # hash one fixed key and could be 0-in-N sampled forever
        key = (args, tuple(sorted(kw.items())), self.calls)
        if self.selected(key):
            self.checked += 1
            o = self.oracle(*args, **kw)
            d = agree_digits(p, o, working_digits=self.digits + 20)
            self.worst = min(self.worst, d)
            if not (d >= self.digits):    # NaN-safe
                raise DualPathDisagreement(
                    f"{self.label}: sampled check at args={args!r} kw={kw!r} "
                    f"agrees to only {d:.1f} digits (bar {self.digits}) — HALT",
                    primary=p, oracle=o, digits=d)
        return p

    def certificate(self):
        return {"label": self.label, "calls": self.calls,
                "checked": self.checked, "worst_agree_digits": self.worst,
                "bar": self.digits, "every": self.every}

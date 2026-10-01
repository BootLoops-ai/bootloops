"""Exact evidence + component-count posterior (front door, small/medium regime).

Delegates to the vendored gated evaluator (jeff/bigg.py) via the integrity-
checked engine loader; g = 1 closed form from the gated identity. For the
scale routes (large k / large g / large N / DP limit) use mixalot.engines
directly — recipes in MANUAL.md."""
from fractions import Fraction
from math import factorial

from .engines import load


def _Z_bigg():
    return load("bigg").Z_bigg


MAX_N_GUARD = 2_000_000     # cost fence: see MANUAL "Costs (measured)"


class NonIntegerCountError(ValueError):
    """Counts must be integers — floats are REJECTED, never
    silently truncated (Z([1.5,2.5],g) must never compute Z([1,2],g))."""


def _validate_counts(U):
    out = []
    for u in U:
        if isinstance(u, bool):
            raise NonIntegerCountError(f"count {u!r} is a bool, not a count")
        if isinstance(u, int):
            v = u
        elif isinstance(u, float):
            if not u.is_integer():
                raise NonIntegerCountError(
                    f"count {u!r} is not an integer; MIXALOT never truncates "
                    "— round your counts explicitly if that is what you mean")
            v = int(u)
        else:
            try:
                from fractions import Fraction as _F
                f = _F(u)
            except Exception:
                raise NonIntegerCountError(f"count {u!r} is not numeric")
            if f.denominator != 1:
                raise NonIntegerCountError(
                    f"count {u!r} is not an integer; MIXALOT never truncates")
            v = int(f)
        if v < 0:
            raise ValueError("counts must be non-negative")
        out.append(v)
    if not out:
        raise ValueError("U must be non-empty")
    return out


def Z(U, g, measure="dirichlet"):
    """Exact evidence Z(U) for a g-component mixture (Fraction).

    U: allocation counts over k states (integers; floats must be integral —
    non-integer counts raise NonIntegerCountError, never truncate).
    measure: "dirichlet" (default, the gated-identity Dir(1) convention) or
    "lebesgue" (divides by (g-1)!*((k-1)!)**g)."""
    U = _validate_counts(U)
    if g < 1:
        raise ValueError("g >= 1")
    if measure not in ("dirichlet", "lebesgue"):
        raise ValueError(
            f"measure must be 'dirichlet' or 'lebesgue', got {measure!r}")
    if len(U) == 1:
        # k = 1: a single state — Z = 1 exactly under both measures' shared
        # normalization of the trivial simplex (m5 shortcut; avoids
        # factorial(N) grinds on huge single-state counts)
        return Fraction(1)
    N = sum(U)
    if N > MAX_N_GUARD:
        raise ValueError(
            f"N={N} exceeds the front-door cost guard ({MAX_N_GUARD:,}); "
            "this regime needs the scale engines (mixalot.engines.bigk2 / "
            "zseries) — see MANUAL 'Costs (measured)'")
    if g == 1:
        k, N = len(U), sum(U)
        z = Fraction(factorial(k - 1), factorial(N + k - 1))
        for u in U:
            z *= factorial(u)
        if measure == "lebesgue":
            z /= factorial(k - 1)
        elif measure != "dirichlet":
            raise ValueError(measure)
        return z
    return _Z_bigg()(U, g, measure)


def bayes_factor(U, g1, g2):
    """Exact BF Z(U,g1)/Z(U,g2) as a Fraction."""
    return Z(U, g1) / Z(U, g2)


def gstar(U, gmax=3, prior=None):
    """Exact posterior over the component count g = 1..gmax (all Fractions;
    sums to 1 by construction). prior: dict g -> positive weight (default
    uniform)."""
    gs = list(range(1, gmax + 1))
    if prior is None:
        prior = {g: Fraction(1, len(gs)) for g in gs}
    else:
        missing = [g for g in gs if g not in prior]
        if missing:
            raise ValueError(f"prior missing weights for g={missing}")
        w = {g: Fraction(prior[g]) for g in gs}
        if any(x <= 0 for x in w.values()):
            raise ValueError(
                "prior weights must be strictly positive (b2: negative/zero "
                f"weights would fabricate pseudo-probabilities); got {prior}")
        tot = sum(w.values())
        prior = {g: w[g] / tot for g in gs}
    ev = {g: Z(U, g) for g in gs}
    joint = {g: prior[g] * ev[g] for g in gs}
    denom = sum(joint.values())
    if denom == 0:
        raise ValueError("zero total evidence (degenerate input)")
    post = {g: joint[g] / denom for g in gs}
    map_g = max(gs, key=lambda g: (post[g], -g))
    return {"posterior": post, "map_g": map_g,
            "p_ge2": sum(post[g] for g in gs if g >= 2),
            "evidence": ev, "prior": prior, "gmax": gmax}


# ---------------- front-door wrappers over vendored engines (input-hardened)

def p_g1(X, gmax=6):
    """Census power probe P(g=1|X) via the vendored engine, with the m1
    guard: the internal grid needs g <= n units (comb(n-1,g-1)>0)."""
    n = len(X)
    if n < 2:
        raise ValueError(f"need at least 2 units, got {n}")
    if n < gmax:
        raise ValueError(
            f"n={n} units supports g only up to {n}; call with gmax<={n} "
            "(the vendored default grid g=1..6 needs n>=6)")
    from .engines import load
    return load("mic_power").p_g1(X)


def lsx_example(name):
    """Runnable route for the LSX flagship examples:
    name in {'swiss','coin10','coin242'}; returns True on PASS."""
    if name not in ("swiss", "coin10", "coin242"):
        raise ValueError("name must be swiss|coin10|coin242")
    from .engines import load
    return bool(load("lsx_direct").run_example(name))


def Zdpm(U, alpha):
    """Exact Dirichlet-process-limit evidence (front door): coerces alpha to
    an exact Fraction so the advertised-exact route always returns the exact
    type, never float (Zdpm_gf((2,1),1.0) must not return a float)."""
    U = _validate_counts(U)
    from .engines import load
    return load("w4_blind_gf").Zdpm_gf(U, Fraction(alpha))

"""baller.frontdoor — the user-facing ball-arithmetic layer (run / solve /
render): BALLER never aliases/shadows Arb's namespace, but it OWNS the
user-facing layer that builds with it. Thin — every ball operation is
python-flint/Arb (arb/acb); nothing here re-implements ball arithmetic.

Entries (all exposed at top level: baller.run / baller.solve / baller.render):

  run(f, x0=None, n=None, dps=50)
      Execute a user expression or iteration in ball arithmetic at working
      precision dps decimal digits (= ceil(dps*log2(10)) bits, no guard bits
      — the working precision IS the user's ask; escalation belongs to
      solve). Expression mode (x0 is None): f is a 0-arg callable or a str
      eval'd with {arb, acb, flint} in scope. Iteration mode: x0 = seed or
      chronological seed tuple (x_0, .., x_{k-1}); f is the step map called
      f(*window) on the last k values OLDEST FIRST, iterated up to index n.
      Muller's recurrence (the reference figure):
          run(lambda xp, x: 111 - 1130/x + 3000/(x*xp), x0=(2, -4), n=100,
              dps=50)
      Returns a RunResult (.values, .value, .blown_at, .render()).

  solve(f, target_digits, ...)
      Adaptive certified driver: run f under ball arithmetic watching the
      radius; on blowup (radius swamping the target, or a non-finite ball)
      escalate the working precision and rerun. ESCALATION POLICY
      (documented, deterministic): start at start_dps (default
      max(target_digits + 10, 30) decimal digits), MULTIPLY BY factor
      (default 2 — precision doubling) after every failed attempt, up to
      max_dps (default 4096); reaching the cap without target_digits
      certified raises SolveRefused — a TYPED refusal NAMING the achieved
      digits, never a bare wrong number. Success returns a SolveResult
      (mid +/- rad with >= target_digits certified by the 1-ulp criterion
      below).

  render(x, digits=None, strict=False)
      FAIL-CLOSED printing: prints ONLY certified digits. A printed string p
      with d significant digits is CERTIFIED iff |p - t| <= 1 ulp of p's
      last printed place for EVERY t in the ball (checked rigorously in arb;
      round-to-nearest display semantics — the reference: a ~300-digit run
      prints x_100 = 6.0000000160995649 against true ...64889..., last
      digit within 1 ulp). Requests beyond the certified count are CLAMPED
      to the certified digits (strict=False) or REFUSED typed
      (strict=True, RenderRefused). Uncertified values (infinite radius,
      non-finite mid, zero-mid fat balls) are returned as a MARKED string
      "UNCERTIFIED[~mid +/- rad]" or refused — never printed bare.

DOCUMENTED DEVIATION from Arb's convention: in the ITERATION driver, when a step's ball goes
non-finite — canonically a divisor ball that CONTAINS 0, where Arb returns
an indeterminate [nan +/- inf] ball — the radius goes to +inf and the
CENTRAL VALUE CONTINUES as the fixed-precision stream: the value at every
subsequent index is arb(mid_i, +inf) where mid_i is the mid-trimmed
iteration at working precision (exactly what plain fixed-precision
arithmetic computes). [anything +/- inf] is a VALID ENCLOSURE of the true
value — this is deliberately NOT Arb's nan/indeterminate choice: the mid
stream keeps telling the fixed-precision story (the Muller float-trap
trajectory to 100) while the infinite radius keeps it honestly uncertified;
render marks it, never prints it bare. One-shot expression mode has no
stream to continue and returns Arb's ball as computed (render still marks a
non-finite result). If the mid stream itself hits an exact zero divisor,
the mid becomes nan and the value is [nan +/- inf] — honest, still marked.

NEGATIVE CONTROL (documented, battery L17): the same Muller recurrence in
plain Python floats converges to EXACTLY 100.0 (true limit 6) — the famous
float trap. That bare wrong number is what fail-closed printing exists to
prevent: the 50-digit ball run marks x_100 UNCERTIFIED (mid ~100, rad inf)
instead, and solve(target 16) returns the certified 6.0000000160995649.

Standalone-loadable by design (no package-relative imports): the battery's
mutation harness (leg L18) gates this file as a {FILE} target.
"""
import math

__all__ = ["run", "solve", "render", "certified_digits", "RunResult",
           "SolveResult", "SolveRefused", "RenderRefused"]

_LOG2_10 = 3.321928094887362


def _bits(dps):
    """Working precision: dps decimal digits -> ceil(dps*log2(10)) bits.
    No guard bits — this entry executes AT the user's stated precision
    (the 50-digit Muller reference radius ~1e-5 is pinned to this map);
    precision escalation is solve's job, never a hidden constant."""
    if not (isinstance(dps, int) and dps >= 1):
        raise ValueError(f"dps must be a positive int, got {dps!r}")
    return int(math.ceil(dps * _LOG2_10))


def _is_acb(v):
    from flint import acb
    return isinstance(v, acb)


def _trim(v):
    """Exact midpoint extraction (the fixed-precision stream's carrier)."""
    from flint import arb, acb
    if isinstance(v, acb):
        return acb(v.real.mid(), v.imag.mid())
    return arb(v.mid())


def _inf_ball(mid_v):
    """The deviation carrier: [mid +/- inf] — a valid enclosure whose
    central value is the fixed-precision stream (NOT Arb's nan choice)."""
    from flint import arb, acb
    if _is_acb(mid_v):
        return acb(arb(mid_v.real.mid(), "inf"), arb(mid_v.imag.mid(), "inf"))
    return arb(mid_v.mid(), "inf")


class RunResult:
    """Result of baller.run. Fields: values (index 0..n in iteration mode,
    length 1 in expression mode), value (final), dps, prec_bits, blown_at
    (first index whose ball leg went non-finite; values from there on are
    the [mid-stream +/- inf] continuation — see the DEVIATION note), mode."""

    def __init__(self, values, dps, prec_bits, blown_at, mode):
        self.values = values
        self.dps = dps
        self.prec_bits = prec_bits
        self.blown_at = blown_at
        self.mode = mode

    @property
    def value(self):
        return self.values[-1]

    def certified_digits(self, index=None):
        v = self.value if index is None else self.values[index]
        return certified_digits(v)

    def render(self, digits=None, index=None, strict=False):
        v = self.value if index is None else self.values[index]
        return render(v, digits=digits, strict=strict)

    def __repr__(self):
        return (f"RunResult(mode={self.mode!r}, dps={self.dps}, "
                f"blown_at={self.blown_at}, value={self.render()})")


def run(f, x0=None, n=None, dps=50):
    """Execute a user expression or iteration in ball arithmetic at working
    precision dps (see module docstring for the full contract, the Muller
    example, and the documented zero-divisor DEVIATION). Thin: every
    operation inside f is python-flint/Arb at ceil(dps*log2(10)) bits."""
    import flint
    from flint import arb, acb, ctx
    prec = _bits(dps)
    old = ctx.prec
    try:
        ctx.prec = prec                     # baller:ctx-ok (restored below)
        if x0 is None:
            if n is not None:
                raise ValueError("n given without seeds x0 — iteration mode "
                                 "needs both")
            if isinstance(f, str):
                val = eval(f, {"__builtins__": {}, "arb": arb, "acb": acb,
                               "flint": flint})
            elif callable(f):
                val = f()
            else:
                raise TypeError("expression mode: f must be a str or a "
                                "0-arg callable")
            return RunResult([val], dps, prec, None, "expression")
        seeds = tuple(x0) if isinstance(x0, (tuple, list)) else (x0,)
        if not callable(f):
            raise TypeError("iteration mode: f must be a callable step map")
        if not (isinstance(n, int) and n >= len(seeds) - 1):
            raise ValueError(f"iteration mode: n must be an int index "
                             f">= {len(seeds) - 1}, got {n!r}")
        cplx = any(isinstance(s, complex) or _is_acb(s) for s in seeds)
        mk = acb if cplx else arb
        nan = (acb(arb("nan"), arb("nan")) if cplx else arb("nan"))
        ball = [mk(s) for s in seeds]       # exact seeds at working prec
        mids = list(ball)                   # the fixed-precision stream
        k = len(seeds)
        values = list(ball)
        blown_at = None
        for i in range(k, n + 1):
            if blown_at is None:
                b = f(*ball[-k:])
                if not isinstance(b, (arb, acb)):
                    b = mk(b)
                if b.is_finite():
                    ball.append(b)
                else:
                    blown_at = i            # divisor ball contained 0 (or
                    #                         overflow): Arb says nan — the
                    #                         DEVIATION takes over below
            try:
                m = f(*mids[-k:])
                if not isinstance(m, (arb, acb)):
                    m = mk(m)
                m = _trim(m)
            except (ZeroDivisionError, ValueError):
                m = nan                     # exact zero divisor in the mid
                #                             stream: honest nan, still inf-rad
            mids.append(m)
            values.append(ball[-1] if blown_at is None else _inf_ball(m))
        return RunResult(values, dps, prec, blown_at, "iteration")
    finally:
        ctx.prec = old                      # baller:ctx-ok (the restore)


# ---------------------------------------------------------------------------
# fail-closed printing


class RenderRefused(ArithmeticError):
    """Typed refusal from render(strict=True): the requested digits are not
    certified. Carries certified (achieved) and requested counts."""

    def __init__(self, msg, certified, requested):
        super().__init__(msg)
        self.certified = certified
        self.requested = requested


def _ulp_exponent(p_str):
    """Exponent of one unit in the LAST printed place of a decimal string
    as printed by arb.str (forms: '6.0056', '-6.0056', '0.000123',
    '1.235e+8', '100.000')."""
    s = p_str.lstrip("+-")
    if "e" in s:
        mant, exp = s.split("e")
        e10 = int(exp)
    else:
        mant, e10 = s, 0
    frac = len(mant.split(".")[1]) if "." in mant else 0
    return e10 - frac


def _certified_str(x, d):
    """The d-significant-digit print of x's midpoint IF certified by the
    1-ulp criterion (|p - t| <= ulp for every t in the ball x, verified as
    a certified arb comparison), else None. The verification precision is
    a deterministic function of d alone (never the ambient ctx), so a
    certified count re-renders identically at any ambient precision."""
    from flint import arb, ctx
    old = ctx.prec
    try:
        ctx.prec = int(d * _LOG2_10) + 64   # baller:ctx-ok (restored below)
        p_str = x.mid().str(d, radius=False, more=True)
        p = arb(p_str)
        ulp = arb(10) ** _ulp_exponent(p_str)
        try:
            certified = bool(abs(p - x) <= ulp)   # True only if provable
        except ValueError:
            certified = False               # undecidable = not certified
        return p_str if certified else None
    finally:
        ctx.prec = old                      # baller:ctx-ok (the restore)


_D_CAP = 100000


def certified_digits(x):
    """Certified significant decimal digits of a ball under the 1-ulp
    display criterion (module docstring). 0 = nothing certifiable (includes
    every infinite-radius / non-finite / zero-mid-fat ball); math.inf for an
    exact ball. acb: the min over the real and imaginary parts (each part
    certified in its own magnitude)."""
    if _is_acb(x):
        return min(certified_digits(x.real), certified_digits(x.imag))
    if not x.is_finite():
        return 0
    r = x.rad()
    if r.is_zero():
        return math.inf
    m = x.mid()
    if m.is_zero():
        return 0                            # no significant digits of 0;
        #                                     absolute bars are RadiusWatch's
    from flint import arb, ctx
    old = ctx.prec
    try:
        ctx.prec = max(old, 96)             # baller:ctx-ok (restored below)
        est = int(float((abs(arb(m)) / arb(r)).log() / arb(10).log())) + 1
    finally:
        ctx.prec = old                      # baller:ctx-ok (the restore)
    for d in range(min(est + 2, _D_CAP), 0, -1):
        if d < est - 60:
            break                           # the estimate is log-accurate;
        #                                     never scan unboundedly
        if _certified_str(x, d) is not None:
            return d
    return 0


def _mark_uncertified(x, strict):
    from flint import ctx
    old = ctx.prec
    try:
        ctx.prec = max(old, 64)             # baller:ctx-ok (restored below)
        if _is_acb(x):
            approx = complex(x.mid())
            rad_s = f"{float(x.real.rad())!r},{float(x.imag.rad())!r}"
        else:
            approx = x.mid().str(8, radius=False, more=True)
            r = x.rad()
            rad_s = (r.str(3, radius=False, more=True) if r.is_finite()
                     else "inf")
    finally:
        ctx.prec = old                      # baller:ctx-ok (the restore)
    msg = f"UNCERTIFIED[~{approx} +/- {rad_s}]"
    if strict:
        raise RenderRefused(f"render: no certified digits — {msg}", 0, None)
    return msg


def render(x, digits=None, strict=False):
    """FAIL-CLOSED printing (contract + the DEVIATION note: module
    docstring). Prints ONLY certified digits; uncertified output is marked
    UNCERTIFIED[...] or refused typed (strict=True) — never printed bare.
    digits=None prints every certified digit (exact balls default to 20)."""
    if isinstance(x, (RunResult, SolveResult)):
        x = x.value if isinstance(x, RunResult) else x.ball
    if digits is not None and not (isinstance(digits, int) and digits >= 1):
        raise ValueError(f"digits must be a positive int, got {digits!r}")
    if not _is_acb(x) and not hasattr(x, "is_finite"):
        raise TypeError(f"render takes arb/acb balls (or Run/SolveResult), "
                        f"got {type(x).__name__} — a bare scalar carries no "
                        f"radius to certify")
    if _is_acb(x):
        if x.imag.is_zero():
            return render(x.real, digits=digits, strict=strict)
        re = render(x.real, digits=digits, strict=strict)
        im = render(x.imag, digits=digits, strict=strict)
        return f"({re} + {im}j)"
    if not x.is_finite():
        return _mark_uncertified(x, strict)
    cert = certified_digits(x)
    if cert == 0:
        return _mark_uncertified(x, strict)
    if cert == math.inf:                    # exact ball: every digit true
        d = digits if digits is not None else 20
        return x.mid().str(d, radius=False, more=True)
    if digits is not None:
        p = _certified_str(x, digits)
        if p is not None:
            return p
        if strict:
            raise RenderRefused(
                f"render: {digits} digits requested, only {cert} certified",
                cert, digits)
        # fail-closed clamp: the certified digits, nothing more
    return _certified_str(x, cert)


# ---------------------------------------------------------------------------
# adaptive certified driver


class SolveRefused(ArithmeticError):
    """Typed refusal from solve: the precision cap was reached with fewer
    than target_digits certified. NAMES the achieved digits (achieved_digits,
    best over all attempts) and carries the attempt log + the best ball —
    never a bare wrong number. The best ball is still a valid enclosure."""

    def __init__(self, msg, achieved_digits, target_digits, attempts, best):
        super().__init__(msg)
        self.achieved_digits = achieved_digits
        self.target_digits = target_digits
        self.attempts = attempts
        self.best = best


class SolveResult:
    """A certified solve: .ball (mid +/- rad, >= target certified), .dps
    (the working precision that succeeded), .certified_digits, .attempts
    [(dps, certified_digits)...], .target_digits. repr = mid +/- rad."""

    def __init__(self, ball, dps, cert, attempts, target_digits):
        self.ball = ball
        self.dps = dps
        self.certified_digits = cert
        self.attempts = attempts
        self.target_digits = target_digits

    @property
    def mid(self):
        return self.ball.mid()

    @property
    def rad(self):
        return self.ball.rad()

    def render(self, digits=None, strict=False):
        return render(self.ball, digits=digits, strict=strict)

    def __repr__(self):
        from flint import ctx
        old = ctx.prec
        try:
            ctx.prec = max(old, 64)         # baller:ctx-ok (restored below)
            rad_s = (self.ball.rad().str(3, radius=False, more=True)
                     if not _is_acb(self.ball) else "(acb)")
        finally:
            ctx.prec = old                  # baller:ctx-ok (the restore)
        return (f"SolveResult({self.render()} +/- {rad_s}, "
                f"certified={self.certified_digits}, dps={self.dps})")


def solve(f, target_digits, x0=None, n=None, index=None,
          start_dps=None, max_dps=4096, factor=2):
    """Adaptive certified driver (full contract: module docstring). Runs f
    under ball arithmetic watching the radius; on blowup — a non-finite
    ball or a radius swamping the target — escalates the working precision
    and reruns. f is EITHER a step map (with x0/n, run through the
    iteration driver; index selects which value to certify, default the
    final one) OR a callable f(dps) returning a ball/RunResult. Returns a
    SolveResult with >= target_digits certified, or raises SolveRefused
    naming the achieved digits. Escalation: start_dps (default
    max(target_digits + 10, 30)) doubled (factor=2) per attempt up to
    max_dps (default 4096)."""
    if not (isinstance(target_digits, int) and target_digits >= 1):
        raise ValueError(f"target_digits must be a positive int, "
                         f"got {target_digits!r}")
    if factor <= 1:
        raise ValueError(f"escalation factor must be > 1, got {factor!r}")
    dps = start_dps if start_dps is not None else max(target_digits + 10, 30)
    if not (isinstance(dps, int) and dps >= 1):
        raise ValueError(f"start_dps must be a positive int, got {dps!r}")
    attempts = []
    best_cert, best_ball, best_dps = -1, None, None
    while True:
        if x0 is not None:
            r = run(f, x0=x0, n=n, dps=dps)
            val = r.values[index] if index is not None else r.value
        else:
            out = f(dps)
            val = out.value if isinstance(out, RunResult) else out
        cert = certified_digits(val)
        cert_n = 0 if cert == 0 else cert
        attempts.append((dps, cert_n))
        if cert_n > best_cert:
            best_cert, best_ball, best_dps = cert_n, val, dps
        if cert >= target_digits:
            return SolveResult(val, dps, cert, attempts, target_digits)
        nxt = int(dps * factor)
        if nxt > max_dps or nxt <= dps:
            raise SolveRefused(
                f"solve: precision cap max_dps={max_dps} reached with "
                f"{best_cert} digits certified (target {target_digits}); "
                f"attempts {attempts} — refusing to return an uncertified "
                f"number (best enclosure carried on the refusal)",
                best_cert, target_digits, attempts, best_ball)
        dps = nxt

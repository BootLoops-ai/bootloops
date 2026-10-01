"""baller.quad — certified quadrature: the POSQ core (tools/posq), aliased in
place by identity. tools/posq remains its own tool with its own battery;
BALLER fronts it and never moves its bytes.
"""
import os

from ._core import alias_registered, TOOLS

_POSQ = os.path.join(TOOLS, "posq")

posq = alias_registered("posq", _POSQ)
kernel_io = alias_registered("kernel_io", _POSQ)

__all__ = ["posq", "kernel_io"]


# --------------------------------------------------------------------------
# acb_integral leg (extends baller.quad).
# Wraps python-flint's rigorous integrator (Arb acb_calc, Petras algorithm)
# behind the house discipline: fail-closed on non-finite/blown balls, the
# analytic-flag LAW surfaced in the signature (the classic sqrt/log misuse
# silently returns a WRONG certified-looking ball — battery leg L15 keeps
# the negative control alive), and a ball -> (mid, rad, digits) reader so
# consumers gate on certified digits, never on the mid alone.
def integral_certified(func, a, b, prec_dps=50, rel_tol_dps=None,
                       max_rad=None, **acb_opts):
    """Certified \\int_a^b func dx as an Arb ball.  func(x, analytic) MUST
    honor the analytic flag for non-meromorphic integrands (sqrt/log branch
    points: return a non-finite ball when analytic=True and analyticity is
    unverified — forward the flag to acb methods that accept it).
    Returns dict(mid_str, rad, certified_digits, ball) — FAIL-CLOSED:
    raises BallBlown if the result is non-finite or rad > max_rad."""
    import flint
    prec_bits = int(prec_dps * 3.3219281) + 16
    old = flint.ctx.prec
    try:
        flint.ctx.prec = prec_bits
        if rel_tol_dps is not None:
            acb_opts.setdefault('rel_tol', flint.arb(10) ** (-rel_tol_dps))
        ball = flint.acb.integral(func, a, b, **acb_opts)
    finally:
        flint.ctx.prec = old
    if not ball.is_finite():
        raise BallBlown('integral_certified: non-finite ball (integrand '
                        'non-analytic on path, pole on path, or analytic '
                        'flag misuse)')
    rad = abs(float(ball.rad()))
    mid = abs(complex(ball.mid()))
    import math
    digits = (math.inf if rad == 0 else
              -math.log10(rad / mid) if mid > 0 else -math.log10(rad))
    if max_rad is not None and rad > max_rad:
        raise BallBlown('integral_certified: rad %.3e > max_rad %.3e '
                        '(halt-not-average law: widen prec or split the '
                        'path, never accept a fat ball)' % (rad, max_rad))
    return {'mid_str': ball.mid().str(int(prec_dps)), 'rad': rad,
            'certified_digits': round(digits, 1), 'ball': ball}


class BallBlown(ArithmeticError):
    """Fail-closed: a certified integral whose ball is non-finite or wider
    than the caller's stated tolerance."""


__all__ = ["posq", "kernel_io", "integral_certified", "BallBlown"]

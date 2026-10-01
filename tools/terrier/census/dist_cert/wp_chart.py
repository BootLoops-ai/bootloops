"""dist_cert.wp_chart -- Weil-Petersson alternative metric, the CROSS-CHECK chart.

Regularization statement (pinned):
  "WP distance = certified line integral of sqrt(g_zzbar) along the straight
   z-segment between ball centers; NO cutoff needed (the conifold singularity
   of g_zzbar is log-divergent but integrable; the distance is finite,
   ~ |dz| sqrt(log 1/|dz|))."
Robustness clause: this chart is the ALTERNATIVE/cross-check only; the
headline chart is the algebraic z coordinate (zchart.py). Conifold-locus
membership is NEVER adjudicated here (it is chart-sensitive).

The g_zzbar model is supplied by the caller as a certified interval extension
(census use: transported periods; tests: SYNTHETIC models only). NO census
data is touched by this module.
"""
from fractions import Fraction

from ball import RIv, sqrt_enclosure, sqrt_iv, ln_enclosure

F = Fraction


class ConstantG:
    """Synthetic model g_zzbar == c (tests / two-chart consistency)."""

    singular_head = None

    def __init__(self, c):
        self.c = F(c)

    def enclosure(self, ta, tb):
        return RIv(self.c, self.c)


class LogG:
    """Synthetic log-divergent model g_zzbar(z(t)) = A*ln(1/(1-t)) + B,
    singular at the t=1 endpoint (conifold-shaped; A>0, A*u+B>=0 on domain).
    enclosure() serves t-panels with tb < 1; the head panel [1-h, 1] is
    integrated by wp_distance via the certified tangent-line bound."""

    def __init__(self, A, B, head=F(1, 64)):
        self.A, self.B = F(A), F(B)
        self.singular_head = (self.A, self.B, F(head))

    def enclosure(self, ta, tb):
        assert tb < 1
        lo = self.B if ta == 0 else self.B + self.A * ln_enclosure(
            1 / (1 - F(ta)))[0]
        hi = self.B + self.A * ln_enclosure(1 / (1 - F(tb)))[1]
        return RIv(max(F(0), lo), hi)


def _head_integral(A, B, h):
    """Certified RIv for I = int_0^h sqrt(A*ln(1/s) + B) ds (A > 0, h < 1).

    Substitution s = e^-u gives I = int_{u0}^inf sqrt(Au+B) e^-u du with
    u0 = ln(1/h) and e^-u0 = h EXACTLY. Upper: tangent line at u0 (sqrt is
    concave) => I <= h*(sqrt(A*u0+B) + A/(2*sqrt(A*u0+B))), evaluated with
    directed enclosures of u0. Lower: integrand of I in s is decreasing,
    so I >= h*sqrt(A*ln(1/h)+B) at the panel edge (u0 lower enclosure).
    This is the no-cutoff integrability bound (~ h*sqrt(ln 1/h)).
    """
    A, B, h = F(A), F(B), F(h)
    assert A > 0 and 0 < h < 1
    u_lo, u_hi = ln_enclosure(1 / h)
    s_hi = sqrt_enclosure(A * u_hi + B)[1]
    s_lo = sqrt_enclosure(max(F(0), A * u_lo + B))[0]
    if s_lo > 0:
        upper = h * (s_hi + A / (2 * s_lo))  # tangent-line (concavity) bound
    else:
        # subadditivity: sqrt(Au+B) <= sqrt(A*u0+B) + sqrt(A(u-u0)) and
        # int_0^inf sqrt(v) e^-v dv = Gamma(3/2) < 1  =>  certified.
        upper = h * (s_hi + sqrt_enclosure(A)[1])
    lower = h * s_lo
    return RIv(lower, upper)


def wp_distance(c1, c2, g_model, n_panels=64, bits=64):
    """Certified RIv for the WP distance between ball CENTERS c1, c2 (CBall),
    per the pinned regularization statement: straight z-segment, certified
    line integral of sqrt(g_zzbar), no cutoff. g_model supplies certified
    interval enclosures of g_zzbar(z(t)), t in [0,1] along the segment
    (plus an optional declared log-singular head at t = 1)."""
    L = (c2 - c1).abs_iv(bits)  # exact-rational |c2-c1| enclosure
    head = getattr(g_model, "singular_head", None)
    t_end = F(1) - (head[2] if head else F(0))
    total = RIv(0, 0)
    for k in range(n_panels):
        ta = t_end * k / n_panels
        tb = t_end * (k + 1) / n_panels
        giv = g_model.enclosure(ta, tb)
        total = total + sqrt_iv(giv, bits) * RIv(tb - ta)
    if head:
        A, B, h = head
        total = total + _head_integral(A, B, h)
    return total * L


def two_chart_consistency(c1, c2, g_model, zdist_iv, n_panels=64, tol=None):
    """Consistency check on synthetic points where BOTH charts are
    defined: the WP enclosure and (scaled) z-chart enclosure must OVERLAP.
    For g == c the exact relation is WP = sqrt(c) * |dz|."""
    w = wp_distance(c1, c2, g_model, n_panels)
    if isinstance(g_model, ConstantG):
        s = sqrt_iv(RIv(g_model.c, g_model.c))
        ref = zdist_iv * s
        return w.overlaps(ref), w, ref
    return w.overlaps(zdist_iv), w, zdist_iv

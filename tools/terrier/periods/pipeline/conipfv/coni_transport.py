#!/usr/bin/env python3
"""
coni_transport.py — (c) certified transport MUM -> conifold-adjacent ->
vacuum locus, reusing the pipe_transport stepped-majorant legs UNCHANGED.

LEG 1 (bulk, z_cf = 0): needs card fields that the shipped coni cards do not carry —
    coni_op     : theta-form bulk operator L_s (same schema as operator_LS;
                  derived by z_cf -> 0 elimination of the card's GKZ boxes,
                  or by annihilating the class-graded restricted series)
    coni_towers : exact log-tower coefficients at the bulk MUM point
    coni_routes : s0/rho0/waypoints (envelope-scan pattern of the DKMM routes)
    s_star      : exact rational landing point
  With those present, this module calls pipe_transport.configure + run_route
  verbatim (the majorant legs are already parametric in the factored
  indicial).  FAIL-CLOSED: every entry point asserts the fields exist.
LEG 2 (transverse jet): d Pi/d z_cf at z_cf = 0 via the pipe_vac pattern
  (connection matrices from the 2-parameter ideal; block-companion local
  legs).  Needs coni_jet_ideal (card input, cross-gated like module_basis).
LEG 3 (conifold-adjacent landing): NOT a series leg.  Ball evaluation of
  the exact Frobenius block at z_cf inside the conifold disc:
    sqrt(pi/2) W = W_bulk + z_cf W1(z_cf) + R,
    |R| <= C2 |z_cf|^2 (1 + |log(-2 pi i z_cf)|)          (*)
  C2 = card input coni_C2 (a certified bound on the z_cf^2 block of the
  2-parameter prepotential; NOT derivable from g1/g2 payload alone).
  (*) is the D2 local form: the z_cf^2-and-higher terms are
  n_cf z^2/(4 pi i) log(-2 pi i z) + analytic; on |z_cf| <= z_max < rho_cf
  the analytic part is bounded by its sup, absorbed into C2 (docstring
  law: C2 must come with its own certificate on the card).
  Ball arithmetic end-to-end (python-flint acb); NEVER reduce through
  float64.
Status: parametric implementation; a card without coni_op makes run()
fail-close (no card in the package carries one).
Deriving coni_op is the coniop/ route (see DESIGN.md SS3).
"""
import sys, os
from fractions import Fraction as Fr

_PIPE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _PIPE)

REQUIRED = ("coni_op", "coni_towers", "coni_routes")


def transport_ready(card):
    """True iff LEG 1 inputs exist on the card (fail-closed router)."""
    return all(card.raw.get(k) for k in REQUIRED) and card.s_star is not None


def run_leg1(card, route, prec):
    """Bulk on-divisor transport: pipe_transport verbatim on coni_op."""
    assert transport_ready(card), \
        "coni_transport LEG1: card lacks coni_op/coni_towers/coni_routes/" \
        "s_star — derive the bulk operator first (DESIGN.md SS3)"
    import pipe_transport as PT
    import pipe_lib as PL
    # configure exactly as run_pipe does for g5, but from the coni_* fields
    raise NotImplementedError(
        "LEG1 wiring activates with the first coni_op-bearing card; the "
        "call is PT.run_route(coni_towers, N_EXT, s0, rho0, waypoints, "
        "targets, prec) — no new transport code (design law).")


def landing_ball(fb, tau_ball, zcf_ball, Wbulk_ball, C2, prec=200):
    """LEG 3: certified sqrt(pi/2) W ball at (tau, z_cf), z_cf in the
    conifold disc.  All inputs acb balls except fb (exact block) and C2
    (exact Fraction upper bound per (*) above)."""
    from flint import acb, arb, ctx
    ctx.prec = prec
    two_pi_i = acb(0, 2) * acb.pi()
    c_log = acb(fb["c_log"].numerator) / fb["c_log"].denominator
    c_tau = acb(fb["c_tau"].numerator) / fb["c_tau"].denominator
    lg = (-two_pi_i * zcf_ball).log()
    W1 = (c_log / two_pi_i) * (lg - 1) + c_tau * tau_ball
    for N, tot, li1_ball in fb.get("li1_balls", []):
        W1 += (acb(tot.numerator) / tot.denominator) * li1_ball / two_pi_i
    W = Wbulk_ball + zcf_ball * W1
    # truncation majorant (*): |R| <= C2 |z|^2 (1 + |log(-2 pi i z)|)
    from pipe_transport import pm_ball    # the suite's +/- ball helper
    zub = abs(zcf_ball).upper()
    lub = abs(lg).upper()
    C2u = arb(C2.numerator) / C2.denominator
    rad = C2u * zub * zub * (arb(1) + lub)
    e = pm_ball(rad)                      # [-rad, rad] enclosure
    return W + acb(e, e)


if __name__ == "__main__":
    from family import load_card
    card = load_card(sys.argv[1])
    print(f"{card.name}: transport_ready = {transport_ready(card)} "
          f"(missing: {[k for k in REQUIRED if not card.raw.get(k)]}"
          f"{' + s_star' if card.s_star is None else ''})")

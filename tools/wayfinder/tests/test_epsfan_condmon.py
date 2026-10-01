"""Condition-monitor tests for vandermonde_laurent_certified:
the 'vand-cond' monitor is LIVE
(cond re-measured at the FINAL post-escalation wp), and the numerically-
singular LU path is NAMED (never an anonymous ZeroDivisionError).
The tests assert the monitor genuinely fires: an escalation
exit that satisfies wp >= dps+guard+ceil(log10 cond_meas)+20 with cond_meas
from the LAST pass has a monitor margin >= 1e-20 BY CONSTRUCTION and an
unreachable raise — a forced saturating cond of 1.0e394 must REFUSE, never
certify silently."""
import os
import sys

import mpmath as mp
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import epsfan  # noqa: E402
from epsfan import vandermonde_laurent_certified, EpsFanCertifyError  # noqa: E402


def f_lau(e):
    x = mp.mpf(e)
    return [3 / x ** 2 - 2 / x + mp.pi + 7 * x]


H1 = ['1/%d' % q for q in [1013, 1117, 1237, 1367, 1511, 1669, 1847, 2039,
                           2251, 2493, 2755, 3037, 3359, 3709, 4099, 4523,
                           4999, 5527, 6101, 6733, 7433, 8209, 9067, 10037]]
H2 = ['1/%d' % q for q in range(11003, 11003 + 36 * 60, 60)]


def test_live_monitor_fires_on_estimator_saturation(monkeypatch):
    # (a) The LIVE monitor genuinely fires — under the cond-estimator
    # saturation measurement model (measured cond ~ 10^(wp+5)), i.e. the
    # the wiring verifier's own defect witness: pre-fix this exact
    # construction CERTIFIED SILENTLY at cond 1.0e394 (verifier_cond_
    # probe2.log B); post-fix it must raise 'vand-cond'.
    #
    # WHY NOT A GENUINE GRID (honesty over theater — measured, e4-condmon-
    # fix/probe_livefire.log{,.v1}): with the wp-escalation in place, a
    # genuine-grid live fire needs the 6-round loop to EXHAUST with the
    # final re-measured cond still violating cond*10^-wp > 10^-(dps+guard),
    # i.e. six successive measurements of the SAME rational node matrix
    # whose ceil(log10 cond) strictly grows, ending >20 orders above the
    # last escalation.  mpmath's exact-rounding arithmetic forbids this:
    # (i) whenever LU completes, the measured cond is order-accurate
    # (differences of rounded near-equal nodes are exact to their leading
    # digits, not noise), so the loop breaks within ~2 rounds with margin
    # >= 1e-20 by its exit condition; (ii) once cluster depth exceeds wp,
    # near-equal nodes round EQUAL and LU collapses to an EXACT zero pivot
    # -> the (now named) singular raise, not a saturating estimate;
    # (iii) the only saturation channel is a ~1-ulp rounding band per node
    # pair, which cannot chain across 6 escalating wp values for a fixed
    # rational grid.  Hence on genuine grids the wp-escalation (+ the named
    # singular raise) SUBSUMES the monitor — as the engine docstring now
    # states — and the live monitor is the fail-closed backstop against
    # estimator measurement pathologies, exercised here via the estimator
    # model itself (measurement stub only; extraction values untouched).
    real_pass = epsfan._vand_pass

    def sat_pass(f_, kmin, kmax, em, spanv, n, wp, eps_nodes=None):
        c, nc, err, _cond = real_pass(f_, kmin, kmax, em, spanv, n, wp,
                                      eps_nodes=eps_nodes)
        with mp.workdps(wp):
            return c, nc, err, mp.mpf(10) ** (wp + 5)

    monkeypatch.setattr(epsfan, '_vand_pass', sat_pass)
    with pytest.raises(EpsFanCertifyError) as exc:
        vandermonde_laurent_certified(f_lau, -2, 1, 30, R0='1/2',
                                      nodes_seed=24, guard=8,
                                      eps_nodes1=H1, eps_nodes2=H2,
                                      max_wp=10 ** 7)
    assert exc.value.kind == 'vand-cond'
    assert exc.value.info.get('singular') is None    # monitor, not LU crash
    assert 'LIVE' in str(exc.value)                  # post-escalation check
    monkeypatch.undo()
    # break path (healthy grid) unchanged: certifies with margin >= 1e-20
    res = vandermonde_laurent_certified(f_lau, -2, 1, 30, R0='1/2',
                                        nodes_seed=24, guard=8,
                                        eps_nodes1=H1, eps_nodes2=H2)
    assert res['diag']['certified'] is True


def test_singular_matrix_raises_named():
    # (b) Ultra-ill-conditioned path: node cluster (rel sep ~1e-18, the
    # wiring-verifier probe config) collapses at working
    # precision.  Pre-fix: anonymous ZeroDivisionError('matrix is
    # numerically singular') — rc!=0 but UNNAMED.  Post-fix: charter-named
    # EpsFanCertifyError('vand-cond') with cond=inf, singular=True and
    # node-separation diagnostics.
    g1 = ['1/%d' % (10 ** 18 + i) for i in range(24)]
    g2 = ['1/%d' % (10 ** 18 + 10 ** 6 + i) for i in range(36)]
    with pytest.raises(EpsFanCertifyError) as exc:
        vandermonde_laurent_certified(f_lau, -2, 1, 30, R0='1/2',
                                      nodes_seed=24, guard=8,
                                      eps_nodes1=g1, eps_nodes2=g2)
    assert exc.value.kind == 'vand-cond'
    assert exc.value.info.get('singular') == 'True'
    assert exc.value.info.get('cond') == 'inf'
    assert 'SINGULAR' in str(exc.value)
    assert 'min_node_sep' in exc.value.info

"""Regression: mplll API calls must NOT depend on ambient mp.dps.

The regression guarded: calling mplll_fit through the API at
ambient mp.dps=15 must not produce junk relations on the rawlll path —
without its own precision block _rawlll_lattice would parse 100d
value-strings and build int(nint(x*10^d)) at 15 digits (the CLI is immune
because mplll_cli sets mp.mp.dps=d+30). Every lattice-build path owns a
mp.workdps(d+guard) block.

Mutation-checked config (scratch-copy revert of the precision block):
rawlll on K=6 functions, 6 points, 4 fit / 2 held-out, planted
I = (314159*F0 - 271828*F1)/100003. Un-patched _rawlll_lattice at ambient
dps=15 returns a height~1e19 junk relation with heldout ~3.7d (the observed
junk-relation mode); easier configs (K=3, or height-7 plants) do NOT discriminate — the
corrupted planted vector stays shortest in the lattice. qrbkz needs
K+1 <= N_fit (LQ inapplicable otherwise) so it runs the K=3 config as an
API-vs-CLI consistency check. In all cases the API call at ambient dps=15
must reproduce the CLI-path (ambient dps = d+30) relation exactly, with
held-out >= 30d.
"""
import mpmath as mp
import pytest

import mplll
from pslq_gate import canonicalize

D = 100        # fit dps (CLI would set ambient to D+30)


def _build(K, npts):
    """(D+30)d value-strings for a K-function basis + planted target."""
    with mp.workdps(D + 60):
        xs = [mp.mpf(3) / 10 + mp.mpf(j) / 7 for j in range(npts)]
        funs = [lambda x: mp.power(mp.pi, x), lambda x: mp.log(2 + x),
                lambda x: x ** 2, lambda x: mp.exp(x),
                lambda x: mp.sqrt(2 + x), lambda x: mp.log(3 + x)]
        F = [[f(x) for x in xs] for f in funs[:K]]
        I = [(314159 * F[0][j] - 271828 * F[1][j]) / 100003 for j in range(npts)]
        s = lambda v: mp.nstr(v, D + 30, strip_zeros=False)
        return [s(v) for v in I], [[s(v) for v in r] for r in F]


# method -> (K, npts, n_fit); rawlll = the mutation-checked ambient-dps config,
# qrbkz = K+1 <= N_fit consistency config (its builder was already wrapped).
CFG = {"rawlll": (6, 6, 4), "qrbkz": (3, 8, 6)}
DATA = {m: _build(k, n) for m, (k, n, _) in CFG.items()}


def _run(method):
    K, npts, nfit = CFG[method]
    I, F = DATA[method]
    names = [f"f{i}" for i in range(K)]
    return mplll.mplll_fit(I, F, names, list(range(nfit)),
                           list(range(nfit, npts)), D,
                           method=method, two_prec=False)


@pytest.mark.parametrize("method", ["rawlll", "qrbkz"])
def test_api_ambient_dps15_matches_cli_path(method):
    K = CFG[method][0]
    want = tuple([100003, -314159, 271828] + [0] * (K - 2))
    old = mp.mp.dps
    try:
        mp.mp.dps = D + 30            # CLI path (mplll_cli sets a.dps+30)
        ref = _run(method)
        mp.mp.dps = 15                # bare API path (ambient dps, the guarded case)
        api = _run(method)
    finally:
        mp.mp.dps = old
    assert ref["relation"] is not None, f"{method}: CLI-path lost planted relation"
    assert canonicalize(ref["relation"]) == want
    assert api["relation"] is not None, \
        f"{method}: API at ambient dps=15 returned NULL (ambient-dps regression)"
    assert canonicalize(api["relation"]) == canonicalize(ref["relation"]), \
        f"{method}: API relation {api['relation']} != CLI relation {ref['relation']}"
    assert api["heldout_min_d"] >= mplll.GATE_D


def test_heldout_digits_ambient_independent():
    I, F = DATA["rawlll"]
    want = [100003, -314159, 271828, 0, 0, 0, 0]
    old = mp.mp.dps
    try:
        mp.mp.dps = 15
        hd15, _ = mplll.heldout_digits(want, I, F, [4, 5], D)
        mp.mp.dps = D + 30
        hdhi, _ = mplll.heldout_digits(want, I, F, [4, 5], D)
    finally:
        mp.mp.dps = old
    assert abs(hd15 - hdhi) < 2 and hd15 >= mplll.GATE_D

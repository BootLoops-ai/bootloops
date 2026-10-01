"""Battery leg: the CLI must accept --method rawbkz (raw lattice + BKZ-beta).

mplll.py's own docstring ships the CLI contract "--method qrbkz|rawlll|rawbkz|cvp"
and GUIDE.md routes users to rawlll/rawbkz when N_points < K+1 (qrbkz's LQ needs
K+1 <= N_fit), but the argparse choices list refused rawbkz — a dead flag exactly
where the guide points. This leg runs the CLI end-to-end in that undersampled
regime (K=6 functions, 4 points, 2 fit / 2 held-out) on a planted relation and
checks the CLI result against the API path. Toy fixture, seconds-scale,
backend-agnostic: with no BKZ backend present, reduce_rows degrades to plain LLL
and records that in the backend string — the planted relation is found either way
at this dimension.
"""
import json
import sys

import mpmath as mp
import pytest

import mplll
import mplll_cli
from pslq_gate import canonicalize

D = 100
K, NPTS = 6, 4          # N_points < K+1: the regime GUIDE routes to rawlll/rawbkz
FIT, HO = [0, 1], [2, 3]
WANT = tuple([100003, -314159, 271828] + [0] * (K - 2))


def _build():
    """(D+30)d value-strings for the K-function basis + planted target."""
    with mp.workdps(D + 60):
        xs = [mp.mpf(3) / 10 + mp.mpf(j) / 7 for j in range(NPTS)]
        funs = [lambda x: mp.power(mp.pi, x), lambda x: mp.log(2 + x),
                lambda x: x ** 2, lambda x: mp.exp(x),
                lambda x: mp.sqrt(2 + x), lambda x: mp.log(3 + x)]
        F = [[f(x) for x in xs] for f in funs[:K]]
        I = [(314159 * F[0][j] - 271828 * F[1][j]) / 100003 for j in range(NPTS)]
        s = lambda v: mp.nstr(v, D + 30, strip_zeros=False)
        return [s(v) for v in I], [[s(v) for v in r] for r in F]


I_STR, F_STR = _build()
NAMES = [f"f{i}" for i in range(K)]


def _run_cli(tmp_path, monkeypatch, capsys, method):
    basis = tmp_path / "basis.json"
    target = tmp_path / "target.json"
    basis.write_text(json.dumps(
        {"names": NAMES, "values": {n: F_STR[i] for i, n in enumerate(NAMES)}}))
    target.write_text(json.dumps(I_STR))
    argv = ["mplll_cli.py", "--basis", str(basis), "--target", str(target),
            "--dps", str(D), "--method", method,
            "--fit-points", ",".join(map(str, FIT)),
            "--held-out", ",".join(map(str, HO))]
    monkeypatch.setattr(sys, "argv", argv)
    old = mp.mp.dps
    try:
        with pytest.raises(SystemExit) as e:
            mplll_cli.main()
    finally:
        mp.mp.dps = old             # main() sets ambient dps; don't leak it
    return e.value.code, json.loads(capsys.readouterr().out)


def test_cli_accepts_rawbkz_in_undersampled_regime(tmp_path, monkeypatch, capsys):
    rc, out = _run_cli(tmp_path, monkeypatch, capsys, "rawbkz")
    assert rc == 0, f"CLI exit {rc}; out: {out}"
    assert out["method"] == "rawbkz" and out["status"] == "HIT"
    assert canonicalize(out["relation"]) == WANT
    assert out["heldout_min_d"] >= mplll.GATE_D
    assert out["two_prec_stable"]


def test_cli_rawbkz_matches_api(tmp_path, monkeypatch, capsys):
    rc, out = _run_cli(tmp_path, monkeypatch, capsys, "rawbkz")
    assert rc == 0, f"CLI exit {rc}; out: {out}"
    api = mplll.mplll_fit(I_STR, F_STR, NAMES, FIT, HO, D,
                          method="rawbkz", two_prec=False)
    assert api["relation"] is not None, "API rawbkz lost the planted relation"
    assert canonicalize(api["relation"]) == canonicalize(out["relation"]) == WANT

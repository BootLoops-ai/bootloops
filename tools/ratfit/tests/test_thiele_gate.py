# Tests for tools/thiele_gate.py
# (a) synthetic rational-function families (several degrees, constant,
#     degenerate-ordering) + gate behavior + CLI;
# (b) positive control: reproduce the reference acheck 88-entry/10032-check
#     result on the archived data (READ-ONLY) and match ACHECK_GATE.json.
import json
import os
import subprocess
import sys
from fractions import Fraction as F

import pytest

import thiele_gate as tg

TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFDIR = os.environ.get('RATFIT_THIELE_REF_DIR', '')  # reference data dir (not shipped; env-pointed): unset/absent -> positive-control tests SKIP


def _sample(f, xs):
    return xs, [f(x) for x in xs]


def _grid(n, den=7):
    return [F(k, den) for k in range(1, n + 1)]


# ------------------------------------------------------------- synthetic fits
def test_constant_family():
    xs, ys = _sample(lambda x: F(3, 4), _grid(12))
    cf = tg.fit(xs, ys)
    assert cf is not None and cf.depth == 1
    assert tg.validate(cf, [(F(99), F(3, 4))])


def test_polynomial_and_rational_degrees():
    funcs = [
        lambda x: 2 * x - F(1, 3),                                # deg (1,0)
        lambda x: (3 * x ** 2 - 1) / (x ** 3 + 7),                # deg (2,3)
        lambda x: (x ** 4 - F(5, 2) * x + 1) / (2 * x ** 2 + 3),  # deg (4,2)
    ]
    for f in funcs:
        xs, ys = _sample(f, _grid(30))
        cf = tg.fit(xs, ys, max_fit=22)     # 8 in-sample held-out points
        assert cf is not None, f
        assert cf.depth < 22, 'no saturation'
        assert tg.validate(cf, _sample_pairs(f, [F(41, 3), F(-8, 5)]))


def _sample_pairs(f, xs):
    return [(x, f(x)) for x in xs]


def test_degenerate_ordering_even_function():
    # f(x)=x^2 on a symmetric grid: interleaved order starts (-3, 3, ...) with
    # equal ys -> zero inverse difference (CFDead) -> the rotation retry must recover.
    f = lambda x: x ** 2
    xs = [F(k) for k in range(-6, 7)]
    ys = [f(x) for x in xs]
    cf = tg.fit(xs, ys, max_fit=9)
    assert cf is not None
    assert tg.validate(cf, _sample_pairs(f, [F(10), F(-11, 2)]))


def test_non_family_returns_none():
    # points not on any low-degree rational function (held-out check fails)
    xs = _grid(20)
    ys = [x ** 2 for x in xs]
    ys[-1] += F(1, 10 ** 12)          # corrupt one point
    assert tg.fit(xs, ys, max_fit=15) is None


def test_gate_catches_corruption():
    f = lambda x: (x ** 2 + 1) / (x + 5)
    xs, ys = _sample(f, _grid(24))
    cf = tg.fit(xs, ys, max_fit=18)
    assert cf is not None
    good = _sample_pairs(f, [F(31, 7), F(-2, 9)])
    bad = [(F(3, 11), f(F(3, 11)) + F(1, 10 ** 30))]
    rep = tg.gate(cf, good + bad)
    assert rep['ok'] == 2 and rep['fail'] == 1
    assert rep['fails'][0]['x'] == str(F(3, 11))


def test_gate_family_and_workers_agree():
    fs = [lambda x: x / (x + 2), lambda x: F(7), lambda x: (1 - x ** 3) / 5]
    xs = _grid(26)
    banked = [(x, [f(x) for f in fs]) for x in xs]
    new = [(x, [f(x) for f in fs]) for x in [F(50, 3), F(-9, 4)]]
    r1 = tg.gate_family(banked, new, max_fit=20)
    r2 = tg.gate_family(banked, new, max_fit=20, workers=3)
    assert r1 == r2
    assert r1['entries_ok'] == 3 and r1['new_ok'] == 6
    assert r1['entries_fail'] == 0 and r1['new_fail'] == 0


def test_duplicate_x_raises():
    with pytest.raises(ValueError):
        tg.fit([F(1), F(1)], [F(2), F(2)])


# --------------------------------------------------------------------- CLI
def test_cli_roundtrip(tmp_path):
    f = lambda x: (2 * x + 1) / (x ** 2 + 3)
    xs = _grid(20)
    spec = {
        'banked': [[str(x), [str(f(x)), str(x + 1)]] for x in xs],
        'new': [[str(F(77, 5)), [str(f(F(77, 5))), str(F(77, 5) + 1)]]],
        'max_fit': 15,
    }
    inp, outp = tmp_path / 'in.json', tmp_path / 'out.json'
    inp.write_text(json.dumps(spec))
    r = subprocess.run([sys.executable, os.path.join(TOOLS, 'thiele_gate.py'),
                        str(inp), '-o', str(outp)], capture_output=True,
                       text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    rep = json.loads(outp.read_text())
    assert rep['entries_ok'] == 2 and rep['new_ok'] == 2
    # corrupt the new point -> exit 1
    spec['new'][0][1][0] = str(f(F(77, 5)) + 1)
    inp.write_text(json.dumps(spec))
    r = subprocess.run([sys.executable, os.path.join(TOOLS, 'thiele_gate.py'),
                        str(inp)], capture_output=True, text=True)
    assert r.returncode == 1


# --------------------------------------- positive control (reference archive)
@pytest.mark.skipif(not os.path.isdir(os.path.join(REFDIR, 'degen/s00_MI20cal')),
                    reason='archived reference data not present')
def test_positive_control_acheck_reproduction():
    """Reproduce acheck.py's reference gate byte-for-byte (READ-ONLY data).

    Node set = DEFULL nodes present when ACHECK_GATE.json was written (165
    nodes were farmed after the gate ran; selected by mtime < gate mtime).
    Must match: entries_ok=88, entries_fail=0, new_ok=10032, new_fail=0 and
    every cf_depth.  ~2-3 min on 16 workers.
    """
    import glob
    D = os.path.join(REFDIR, 'degen/s00_MI20cal')
    ref = json.load(open(os.path.join(REFDIR, 'ACHECK_GATE.json')))
    gate_t = os.path.getmtime(os.path.join(REFDIR, 'ACHECK_GATE.json'))
    stage1 = {os.path.basename(f)[4:-5]
              for f in glob.glob(D + '/DE_u*.json')}
    banked, new = [], []
    for fn in sorted(glob.glob(D + '/DEFULL_u*.json')):
        if os.path.getmtime(fn) >= gate_t:
            continue
        tag = os.path.basename(fn)[8:-5]
        n, d = tag.split('_')
        u = F(int(n), int(d))
        J = json.load(open(fn))
        A = [[F(x) for x in row] for row in J['A_full_exact']]
        (banked if (tag in stage1 and J.get('banked_A_match'))
         else new).append((u, A))
    banked.sort()
    new.sort()
    assert len(banked) == 175 and len(new) == 114
    nr, nc = 8, len(banked[0][1][0])
    labels = [f'{i},{j}' for i in range(nr) for j in range(nc)]
    flat_b = [(u, [A[i][j] for i in range(nr) for j in range(nc)])
              for u, A in banked]
    flat_n = [(u, [A[i][j] for i in range(nr) for j in range(nc)])
              for u, A in new]
    rep = tg.gate_family(flat_b, flat_n, max_fit=140, labels=labels,
                         workers=16)
    assert rep['entries_ok'] == ref['entries_ok'] == 88
    assert rep['entries_fail'] == ref['entries_fail'] == 0
    assert rep['new_ok'] == ref['new_ok'] == 10032
    assert rep['new_fail'] == ref['new_fail'] == 0
    assert rep['cf_depths'] == ref['cf_depths']

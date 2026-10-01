"""The CY4 (1,1,1,1,1,25) c_{3/2} fixture leg of the coalescer closure legs (gate_compare.py + pslq_calpha.py
through lockpick.pslq_gate) on the CY4 inputs of record, beside the CY3 leg in test_pslq_gate_route.py.

Fixtures (tests/fixtures/, sha256-pinned in PINS.json): the Route A dps-80 output ROUTEA_CY4_dps80.json, the two
Route B outputs at distinct s_dec (ROUTEB_CY4_1280_sdec14.json at 1280 bits, s_dec 1/4; ROUTEB_CY4_896_sdec740.json
at 896 bits, s_dec 7/40), the K3 (1,1,1,9) control files shared with the CY3 leg, and the records these tools
produced on them: GATE_CY4.json (gate_compare --family CY4 --closed auto), gate_cy4.log (its stdout after the
mask below) and PSLQ_CY4.json (pslq_calpha --family CY4 --from-gate GATE_CY4.json --ring extended).  Every
expected vector and digit count is read from those records (canonicalized) AND asserted as the record's literals
beside them (248 / 261 / 870, 248.22, 70.53, 70.26, 257.93, 74.8 / 223.28, 195, the fit pairs (178,208) / (125,155), the
selftest count 25): a literal that drifts from its pinned fixture fails the leg by name.

The targets rule this leg pins: GATE.json "targets" carries the best Route B midpoint trimmed to
int(min(floor(acc_bits*log10 2), B s_dec-independence)) digits when two or more distinct s_dec are present
(261 ball digits and 248.22 d of s_dec-independence give 248 at CY4), and to the ball digits alone with a printed
CAVEAT line when only one s_dec is present.
"""
import hashlib, json, math, os, re, subprocess, sys

import mpmath as mp
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
FIX = os.path.join(HERE, 'fixtures')
PY = sys.executable
sys.path.insert(0, PKG)
import calpha_rings as R  # noqa: E402

MASK = ('TARGET ', 'GATE_VERDICT', 'PRODUCER')
CY4_FIXTURES = ['ROUTEA_CY4_dps80.json', 'ROUTEB_CY4_1280_sdec14.json', 'ROUTEB_CY4_896_sdec740.json',
                'GATE_CY4.json', 'gate_cy4.log', 'PSLQ_CY4.json']
GATE_VERDICT_OF_RECORD = 'GATE_VERDICT family=CY4 alphas=3/2 cross_route_min_d=70.53 bar_d=30.0 closed=yes k3_control=yes PASS'


def canon(vec):
    g = 0
    for v in vec: g = math.gcd(g, abs(int(v)))
    out = [int(v) // g for v in vec]
    for v in out:
        if v:
            if v < 0: out = [-q for q in out]
            break
    return tuple(out)


def fx(name):
    return os.path.join(FIX, name)


def run(script, args, cwd):
    p = subprocess.run([PY, os.path.join(PKG, script)] + args, cwd=cwd, capture_output=True, text=True,
                       env={'PATH': os.environ.get('PATH', '/usr/bin:/bin'), 'PYTHONDONTWRITEBYTECODE': '1'})
    return p.returncode, p.stdout + p.stderr


def masked(text):
    m = '\n'.join(l for l in text.split('\n') if not l.startswith(MASK))
    return re.sub(r'wrote \S*GATE\.json', 'wrote GATE.json', m)


def ndigits(s):
    return len(s.lstrip('-').replace('.', '').lstrip('0'))


def agree_d(x, y, dps=320):
    with mp.workdps(dps):
        x, y = mp.mpf(x), mp.mpf(y)
        if x == y: return float(dps)
        return float(-mp.log10(abs(x - y) / abs(y)))


GATE_CY4 = json.load(open(fx('GATE_CY4.json')))
PSLQ_CY4 = json.load(open(fx('PSLQ_CY4.json')))
SCANS = {(x['form'], x['alpha']): x for x in PSLQ_CY4['scans']}
GATE_ARGS = ['--family', 'CY4', '--routeA', fx('ROUTEA_CY4_dps80.json'),
             '--k3-routeA', fx('ROUTEA_K3.json'), '--k3-routeB', fx('ROUTEB_K3_ctrl.json')]
TWO_SDEC = f"{fx('ROUTEB_CY4_1280_sdec14.json')},{fx('ROUTEB_CY4_896_sdec740.json')}"


# ------------------------------------------------------------------ integrity
def test_cy4_fixtures_present_and_pinned():
    pins = json.load(open(fx('PINS.json')))
    for f in CY4_FIXTURES:
        assert f in pins, f
        assert hashlib.sha256(open(fx(f), 'rb').read()).hexdigest() == pins[f], f


def test_pf_cy4_data_bank_beside_the_k3_and_cy3_banks():
    p = os.path.join(PKG, 'pf_cy4_data.jl')
    assert os.path.exists(p) and os.path.exists(os.path.join(PKG, 'pf_k3_data.jl')) and os.path.exists(os.path.join(PKG, 'pf_cy3_data.jl'))
    head = open(p).read(4000)
    assert head.startswith('# auto-generated: msq=[1, 1, 1, 1, 1, 25] order=8 degz=14')
    assert 'const ORDER = 8' in head and 'const MSQ = [1, 1, 1, 1, 1, 25]' in head and 'const FRAC_EXPS = [3//2]' in open(p).read()


# ------------------------------------------------------------------ the ring table
def test_cy4_family_table_carries_the_recorded_structure():
    spec = R.family_spec('CY4')
    assert spec['pi_power'] == '2'                                   # the additive-form power the record found
    cf = R.closed_form(spec, '3/2')
    assert cf is not None
    assert cf.to_dict() == {'sign': -1, 'rational': '1/40', 'pi': '-2', 'radicals': {'5': '1/2'}, 'gamma': {},
                            'display': '-sqrt(5)/(40pi^2)'} and cf.expr == '-sqrt(5)/(40*pi^2)'
    names = R.ring_names(spec, 'log', 'extended')
    rec = SCANS[('log', '3/2')]
    assert names == rec['names']                                     # the same basket the record scanned
    assert cf.log_vector(names) == canon(rec['HIT']['vector'])       # the vector DERIVED from the table == the record's HIT
    with mp.workdps(320):
        assert agree_d(cf.value(), GATE_CY4['targets']['3/2']['value']) >= GATE_CY4['targets']['3/2']['digits'] - 1
        assert agree_d(R.parse_expr(cf.expr), cf.value()) >= 300


def test_rings_selftest_count_moves_with_the_table():
    rc, text = run('calpha_rings.py', [], HERE)
    n_closed = sum(len(R.family_spec(f)['closed']) for f in R.families())
    m = re.search(r'^CALPHA_RINGS SELFTEST PASS \((\d+)/(\d+), 0 skipped\)$', text, re.M)
    assert rc == 0 and m, text
    assert int(m.group(1)) == int(m.group(2)) == 13 + 3 * n_closed == 25   # 13 ring facts + 3 checks per recorded closed form (4 forms)


# ------------------------------------------------------------------ the CY4 fixture: gate
def test_gate_compare_cy4_reproduces_the_record(tmp_path):
    out = tmp_path / 'GATE.json'
    rc, text = run('gate_compare.py', GATE_ARGS + ['--routeB', TWO_SDEC, '--out', str(out)], tmp_path)
    assert rc == 0, text
    assert masked(text) == open(fx('gate_cy4.log')).read()           # byte-identical to the log of record after the mask
    assert GATE_VERDICT_OF_RECORD in text
    g = json.load(open(out))
    for k in GATE_CY4:
        if k != 'producer':
            assert g[k] == GATE_CY4[k], k                             # every record key, value for value
    t = g['targets']['3/2']
    assert t['digits'] == 248 and t['ball_digits'] == 261 and t['acc_bits'] == 870 and t['B_sdec_ind_d'] == 248.22
    assert t['digits'] == int(min(t['ball_digits'], t['B_sdec_ind_d'])) and ndigits(t['value']) == 248
    assert t['rule'].startswith('int(min(') and 'CAVEAT' not in t['rule'] and t['source'] == 'ROUTEB_CY4_1280_sdec14.json'
    per = g['c_3/2']['per_sample_vs_closed']
    assert per['cross_route_gate_d'] == 70.53 and per['A_sdec_ind_d'] == 70.26 and per['B_sdec_ind_d'] == 248.22
    assert per['B:B@1//4'] == 248.22 and per['B:B@7//40'] == 257.93
    assert g['K3_control'] == {'A:K3_control@0.15': 74.8, 'A:K3_control@0.25': 74.8, 'B:B@1//4': 223.28}
    assert g['gate'] == {'bar_d': 30.0, 'cross_route_gate_d': {'3/2': 70.53}, 'verdict': 'PASS', 'fails': []}
    assert re.search(r'^TARGET c_\{3/2\}: 248 d \[ROUTEB_CY4_1280_sdec14\.json, acc_bits 870\] rule: int\(min\(', text, re.M)


def test_gate_compare_cy4_single_sdec_prints_the_caveat(tmp_path):
    out = tmp_path / 'GATE.json'
    rc, text = run('gate_compare.py', GATE_ARGS + ['--routeB', fx('ROUTEB_CY4_1280_sdec14.json'), '--out', str(out)], tmp_path)
    assert rc == 0, text                                              # the cross-route gate stands; only the target moves
    g = json.load(open(out))
    t = g['targets']['3/2']
    assert 'B_sdec_ind_d' not in g['c_3/2']['per_sample_vs_closed'] and t['B_sdec_ind_d'] is None
    assert t['digits'] == 261 == t['ball_digits'] and ndigits(t['value']) == 261
    assert 'CAVEAT single s_dec: the ball alone bounds this string; no s_dec-independence member' in t['rule']
    assert re.search(r'^TARGET c_\{3/2\}: 261 d \[ROUTEB_CY4_1280_sdec14\.json, acc_bits 870\] rule: ball digits .* CAVEAT single s_dec', text, re.M)
    assert agree_d(t['value'], GATE_CY4['targets']['3/2']['value']) >= 247   # the same midpoint, longer


# ------------------------------------------------------------------ the CY4 fixture: PSLQ from the gate
@pytest.fixture(scope='module')
def cy4_pslq(tmp_path_factory):
    d = tmp_path_factory.mktemp('pslq_cy4')
    out = d / 'PSLQ.json'
    rc, text = run('pslq_calpha.py', ['--family', 'CY4', '--from-gate', fx('GATE_CY4.json'), '--k3-routeB', fx('ROUTEB_K3_ctrl.json'),
                                      '--ring', 'extended', '--form', 'log,additive', '--out', str(out)], d)
    return rc, text, json.load(open(out))


def test_pslq_cy4_from_gate_refinds_the_record(cy4_pslq):
    rc, text, j = cy4_pslq
    assert rc == 0, text
    assert 'PIN OK' in text and j['status'] == 'CLOSED' and j['pi_power'] == '2'
    assert j['targets']['3/2']['digits'] == 248 and j['targets']['K3_control']['digits'] == 195
    for key, rec in SCANS.items():
        s = {(x['form'], x['alpha']): x for x in j['scans']}[key]
        assert s['names'] == rec['names']
        assert tuple(s['HIT']['vector']) == canon(rec['HIT']['vector']) and s['HIT']['height'] == rec['HIT']['height']
        assert s['verdict'] == 'CLOSED' and s['expectation']['match'] is True and s['HIT']['reverify']['ok']
        assert s['dps_plan']['fit_pair'] == rec['dps_plan']['fit_pair'] and s['dps_plan']['held_out_digits'] >= 35
    assert SCANS[('log', '3/2')]['dps_plan']['fit_pair'] == [178, 208] and SCANS[('log', 'K3_control')]['dps_plan']['fit_pair'] == [125, 155]
    assert 'x pi^2' in {x['alpha']: x for x in j['scans'] if x['form'] == 'additive'}['3/2']['expectation']['note']
    assert j['verdicts'] == {'log:3/2': 'CLOSED', 'log:K3_control': 'CLOSED', 'additive:3/2': 'CLOSED'}


def test_pslq_cy4_controls_inside_the_basket(cy4_pslq):
    rc, text, j = cy4_pslq
    for form, ctls in j['controls'].items():
        for pair, c in ctls.items():
            assert c['ok'] and c['negative']['ok'] and c['negative']['got'] is None
            assert any(d['control'].startswith('planted') for d in c['detail'])


# ------------------------------------------------------------------ controls (each fails BY NAME)
def wrong_closed(tmp_path):
    p = tmp_path / 'wrong_closed.json'
    json.dump({"3/2": {"sign": -1, "rational": "1/36", "pi": "-2", "radicals": {"5": "1/2"}, "gamma": {},
                       "display": "WRONG: -sqrt(5)/(36pi^2)"}}, open(p, 'w'))
    return str(p)


def test_control_wrong_closed_form_fails_the_gate_by_name(tmp_path):
    out = tmp_path / 'GATE.json'
    rc, text = run('gate_compare.py', GATE_ARGS + ['--routeB', TWO_SDEC, '--closed', wrong_closed(tmp_path), '--out', str(out)], tmp_path)
    assert rc == 1
    g = json.load(open(out))
    assert g['gate']['verdict'] == 'FAIL' and g['gate']['cross_route_gate_d'] == {'3/2': 70.53}   # the routes still agree; the form is wrong
    assert re.search(r'c_\{3/2\} B:B@1//4 vs closed \d+\.\d+ d < bar 30\.0', text) and 'closed=yes k3_control=yes FAIL' in text
    assert all(v < 30 for k, v in g['c_3/2']['per_sample_vs_closed'].items() if ' vs ' not in k and k.split(':')[0] in ('A', 'B'))


def test_control_wrong_closed_form_is_a_named_pslq_mismatch(tmp_path):
    wrong = wrong_closed(tmp_path)
    rc, text = run('pslq_calpha.py', ['--family', 'CY4', '--from-gate', fx('GATE_CY4.json'), '--ring', 'extended', '--form', 'log',
                                      '--closed', wrong, '--out', str(tmp_path / 'c.json')], tmp_path)
    names = R.ring_names(R.family_spec('CY4'), 'log', 'extended')
    expected = R.ClosedForm.from_dict(json.load(open(wrong))['3/2']).log_vector(names)
    assert rc == 1 and f'EXPECTATION MISMATCH: expected {expected}' in text and f'got {canon(SCANS[("log", "3/2")]["HIT"]["vector"])}' in text


def test_control_fixture_byte_flip_is_caught_by_the_pin_and_by_the_gate(tmp_path):
    src = open(fx('ROUTEB_CY4_1280_sdec14.json')).read()
    m = re.search(r'"c_re_mid":\s*"\[-0\.00566402635462', src)
    assert m
    i = m.end() - 1
    flipped = src[:i] + ('3' if src[i] != '3' else '4') + src[i + 1:]
    assert flipped != src
    p = tmp_path / 'ROUTEB_CY4_1280_sdec14.json'
    p.write_text(flipped)
    pins = json.load(open(fx('PINS.json')))
    assert hashlib.sha256(flipped.encode()).hexdigest() != pins['ROUTEB_CY4_1280_sdec14.json']   # test_fixture_pins would fail on it
    out = tmp_path / 'GATE.json'
    rc, text = run('gate_compare.py', GATE_ARGS + ['--routeB', f"{p},{fx('ROUTEB_CY4_896_sdec740.json')}", '--out', str(out)], tmp_path)
    assert rc == 1 and 'B:B@1//4 vs closed' in text and '< bar 30.0' in text and 'FAIL' in text
    g = json.load(open(out))
    assert g['c_3/2']['per_sample_vs_closed']['B:B@1//4'] < 30 and g['c_3/2']['per_sample_vs_closed']['B_sdec_ind_d'] < 30
    assert g['c_3/2']['per_sample_vs_closed']['B:B@7//40'] == 257.93            # the untouched sample stands

"""The coalescer closure legs (gate_compare.py + pslq_calpha.py through lockpick.pslq_gate) on the CY3
(1,1,1,1,16) fixture of record and the named controls.

Fixtures (tests/fixtures/, sha256-pinned in PINS.json): the Route A / Route B outputs of the CY3 case and
the K3 (1,1,1,9) control, the PSLQ record PSLQ_PRODUCTION.json, the gate log of record gate.log; the gate
record GATE.json is the one shipped beside the package (tools/coalescer/GATE.json).  Every expected vector
is read from PSLQ_PRODUCTION.json (canonicalized) — nothing is typed here.  The full record dps pair
(120, 200) runs in about a second, so no reduced tier is needed.
"""
import hashlib, json, math, os, re, subprocess, sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
FIX = os.path.join(HERE, 'fixtures')
PY = sys.executable
sys.path.insert(0, PKG)
import calpha_rings as R  # noqa: E402
import pslq_calpha as P   # noqa: E402


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


RECORD = json.load(open(fx('PSLQ_PRODUCTION.json')))
GATE_RECORD = json.load(open(os.path.join(PKG, 'GATE.json')))


# ------------------------------------------------------------------ integrity
def test_fixture_pins():
    pins = json.load(open(fx('PINS.json')))
    for f, h in pins.items():
        assert hashlib.sha256(open(fx(f), 'rb').read()).hexdigest() == h, f
    assert set(pins) == {f for f in os.listdir(FIX) if f != 'PINS.json'}


def test_engine_protocol_pin_matches_the_package_engine():
    """Drift guard: the protocol units of tools/lockpick/pslq_gate.py hash to the pin pslq_calpha carries.
    A deliberate engine change re-pins here (captures the change by name), never silently."""
    eng = os.path.join(PKG, '..', 'lockpick', 'pslq_gate.py')
    assert P.protocol_sha(eng) == P.PSLQ_GATE_PROTOCOL_SHA256


def test_rings_selftest():
    assert R.selftest(dps=100, verbose=False) == 0


# ------------------------------------------------------------------ the CY3 fixture: gate
def test_gate_compare_cy3_reproduces_the_record(tmp_path):
    out = tmp_path / 'GATE.json'
    rc, text = run('gate_compare.py', ['--family', 'CY3', '--routeA', fx('ROUTEA_CY3.json'),
                                       '--routeB', f"{fx('ROUTEB_CY3_sdec025.json')},{fx('ROUTEB_CY3_sdec0175.json')}",
                                       '--k3-routeA', fx('ROUTEA_K3.json'), '--k3-routeB', fx('ROUTEB_K3_ctrl.json'),
                                       '--out', str(out)], tmp_path)
    assert rc == 0, text
    masked = '\n'.join(l for l in text.split('\n') if not l.startswith(('TARGET ', 'GATE_VERDICT', 'PRODUCER')))
    masked = re.sub(r'wrote \S*GATE\.json', 'wrote GATE.json', masked)
    assert masked == open(fx('gate.log')).read()      # byte-identical to the log of record after the mask
    g = json.load(open(out))
    for k in GATE_RECORD:
        assert g[k] == GATE_RECORD[k], k                  # every record key, value for value
    assert g['gate']['verdict'] == 'PASS' and g['gate']['cross_route_gate_d'] == {'5/4': 72.53, '7/4': 71.02} or \
        g['gate']['cross_route_gate_d'] == {al: GATE_RECORD[f'c_{al}']['per_sample_vs_closed']['cross_route_gate_d'] for al in ('5/4', '7/4')}
    assert 'GATE_VERDICT family=CY3 alphas=5/4,7/4 cross_route_min_d=71.02 bar_d=30.0 closed=yes k3_control=yes PASS' in text
    assert g['targets']['5/4']['digits'] == 198 and g['targets']['7/4']['digits'] == 197   # 660 / 657 acc bits -> floor(bits*log10 2)


def test_gate_compare_refusals(tmp_path):
    wrong = tmp_path / 'outside.json'
    json.dump({"5/4": {"sign": -1, "rational": "1", "pi": "-1/2", "radicals": {"2": "-1/2"}, "gamma": {"1/5": -2}}}, open(wrong, 'w'))
    rc, text = run('gate_compare.py', ['--family', 'CY3', '--closed', str(wrong), '--routeA', fx('ROUTEA_CY3.json'),
                                       '--routeB', fx('ROUTEB_CY3_sdec025.json'), '--out', str(tmp_path / 'g.json')], tmp_path)
    assert rc == 3 and 'REFUSED: NotRingMember: closed form needs logG(1/5)' in text
    rc, text = run('gate_compare.py', ['--family', 'CY3', '--routeA', fx('ROUTEA_CY3.json'),
                                       '--routeB', fx('MISSING.json'), '--out', str(tmp_path / 'g.json')], tmp_path)
    assert rc == 3 and 'REFUSED: named input missing' in text
    rc, text = run('gate_compare.py', ['--family', 'CY3', '--routeA', fx('ROUTEA_CY3.json'),
                                       '--routeB', fx('ROUTEB_CY3_sdec025.json'), '--bar', '80', '--out', str(tmp_path / 'g.json')], tmp_path)
    assert rc == 1 and 'cross-route gate 72.53 d < bar 80.0' in text


# ------------------------------------------------------------------ the CY3 fixture: PSLQ at the record's dps pair
@pytest.fixture(scope='module')
def cy3_pslq(tmp_path_factory):
    d = tmp_path_factory.mktemp('pslq')
    out = d / 'PSLQ.json'
    rc, text = run('pslq_calpha.py', ['--family', 'CY3', '--from-routeB', fx('ROUTEB_CY3_hi.json'),
                                      '--k3-routeB', fx('ROUTEB_K3_ctrl.json'), '--dps-pair', '120,200',
                                      '--maxcoeff', '1e6', '--maxcoeff-additive', '1e8', '--out', str(out)], d)
    return rc, text, json.load(open(out))


def test_pslq_engine_bound_by_pin(cy3_pslq):
    rc, text, j = cy3_pslq
    assert rc == 0, text
    e = j['producer']['engine']
    assert e['protocol_sha256'] == e['pinned_protocol_sha256'] == P.PSLQ_GATE_PROTOCOL_SHA256
    assert e['label'] != 'UNLISTED' and 'PIN OK' in text


@pytest.mark.parametrize('form,alpha,record_key', [('log', '5/4', 'log_c5/4'), ('log', '7/4', 'log_c7/4'),
                                                   ('log', 'K3_control', 'K3_log_c32'),
                                                   ('additive', '5/4', 'add_c5/4_piexp'), ('additive', '7/4', 'add_c7/4_piexp')])
def test_pslq_cy3_vectors_match_the_record(cy3_pslq, form, alpha, record_key):
    rc, text, j = cy3_pslq
    s = {(x['form'], x['alpha']): x for x in j['scans']}[(form, alpha)]
    r = RECORD['layers'][record_key]
    assert ['T'] + s['names'] == r['basis']                         # the same basket names
    assert tuple(s['HIT']['vector']) == canon(r['relation'])         # the same canonical vector
    assert s['verdict'] == 'CLOSED' and s['expectation']['match'] is True and s['HIT']['reverify']['ok']
    assert s['dps_plan']['held_out_digits'] >= 30


def test_pslq_cy3_targets_are_the_record_strings(cy3_pslq):
    rc, text, j = cy3_pslq
    assert j['targets']['5/4']['head40'] == RECORD['c54'][:42] and j['targets']['7/4']['head40'] == RECORD['c74'][:42]
    assert j['targets']['5/4']['digits'] == 314 and j['targets']['7/4']['digits'] == 313   # 1044 / 1041 acc bits
    assert {x['alpha']: x['dps_plan']['fit_pair'] for x in j['scans'] if x['alpha'] != 'K3_control'} == {'5/4': [120, 200], '7/4': [120, 200]}
    assert RECORD['layers']['log_c5/4']['dps_tested'] == [120, 200]


def test_pslq_controls_inside_the_basket(cy3_pslq):
    rc, text, j = cy3_pslq
    for form, ctls in j['controls'].items():
        for pair, c in ctls.items():
            assert c['ok'] and c['negative']['ok'] and c['negative']['got'] is None
            assert any(d['control'].startswith('planted') for d in c['detail'])


def test_pslq_k3_family_additive_height_36(tmp_path):
    out = tmp_path / 'K3.json'
    rc, text = run('pslq_calpha.py', ['--family', 'K3', '--from-routeB', fx('ROUTEB_K3_ctrl.json'), '--out', str(out)], tmp_path)
    assert rc == 0, text
    j = json.load(open(out))
    s = {(x['form'], x['alpha']): x for x in j['scans']}
    assert s[('additive', '3/2')]['HIT']['vector_named'] == {'T': 36, '1': 0, 'sqrt3': 1}
    assert tuple(s[('log', '3/2')]['HIT']['vector']) == canon(RECORD['layers']['K3_log_c32']['relation'])[:7]


# ------------------------------------------------------------------ the four controls (each fails BY NAME)
def test_control_planted_wrong_constant_gives_null_never_a_spurious_relation(tmp_path):
    out = tmp_path / 'c1.json'
    rc, text = run('pslq_calpha.py', ['--family', 'CY3', '--from-routeB', fx('ROUTEB_CY3_hi.json'), '--dps-pair', '120,200',
                                      '--form', 'log', '--planted-wrong', 'logG(1/4)=log(gamma(1/3))', '--out', str(out)], tmp_path)
    assert rc == 1
    assert 'NULL-with-floor' in text and 'EXPECTATION MISMATCH: expected (2, 0, 1, 1, 0, 0, 0, 4)' in text and 'got None' in text
    j = json.load(open(out))
    assert all('HIT' not in s for s in j['scans']) and all(c['ok'] for f in j['controls'].values() for c in f.values())


def test_control_degenerate_basket_raises_DegeneratePoolError_by_name(tmp_path):
    out = tmp_path / 'c3.json'
    rc, text = run('pslq_calpha.py', ['--family', 'CY3', '--from-routeB', fx('ROUTEB_CY3_hi.json'), '--dps-pair', '120,200',
                                      '--form', 'log', '--add-member', 'logG(3/4)=log(gamma(3/4))', '--out', str(out)], tmp_path)
    assert rc == 3 and 'REFUSED: DegeneratePoolError: degenerate pool' in text and 'NULL-with-floor' not in text


def test_control_short_target_raises_DigitsError_by_name(tmp_path):
    rc, text = run('pslq_calpha.py', ['--family', 'CY3', '--c', '5/4=' + RECORD['c54'][:62], '--form', 'log',
                                      '--out', str(tmp_path / 'c4.json')], tmp_path)
    assert rc == 4 and re.search(r'^DigitsError: target 5/4 carries only \d+ digits; refusing \(need >= 150\)', text, re.M)


def test_control_capacity_overflow_raises_CapacityError_by_name(tmp_path):
    rc, text = run('pslq_calpha.py', ['--family', 'CY3', '--from-routeB', fx('ROUTEB_CY3_hi.json'), '--dps-pair', '30,200',
                                      '--form', 'log', '--out', str(tmp_path / 'c5.json')], tmp_path)
    assert rc == 4 and re.search(r'^CapacityError: capacity \d+d > available 30d', text, re.M)


def test_control_unlisted_engine_is_refused_by_default_and_the_pin_is_the_second_gate(tmp_path):
    """An engine whose file sha is not in KNOWN_ENGINE_FILES is refused by default (rc 3, by name); under
    --allow-unlisted-engine the protocol pin still gates: a protocol-unit tamper is PIN MISMATCH rc 3, a
    docstring-only edit runs with the loud UNLISTED warning line.  The retired --strict-engine is a usage error."""
    src = open(os.path.join(PKG, '..', 'lockpick', 'pslq_gate.py')).read()
    bad = src.replace('            if v < 0: out = [-q for q in out]\n', '            if v < 0: pass\n')
    assert bad != src
    eng = tmp_path / 'pslq_gate_tampered.py'
    eng.write_text(bad)
    base = ['--family', 'CY3', '--from-routeB', fx('ROUTEB_CY3_hi.json'), '--dps-pair', '120,200', '--form', 'log']
    rc, text = run('pslq_calpha.py', base + ['--engine', str(eng), '--out', str(tmp_path / 'c6.json')], tmp_path)
    assert rc == 3 and 'REFUSED: engine file sha256' in text and 'not in KNOWN_ENGINE_FILES; protocol units DIFFER FROM the pin' in text
    assert 'UNLISTED FILE REFUSED' in text and 'HIT at' not in text
    rc, text = run('pslq_calpha.py', base + ['--engine', str(eng), '--allow-unlisted-engine', '--out', str(tmp_path / 'c6a.json')], tmp_path)
    assert rc == 3 and 'REFUSED: pslq_gate protocol pin mismatch' in text and 'PIN MISMATCH' in text and 'HIT at' not in text
    doc = src.replace('"""pslq_gate.py — shared PSLQ closure harness', '"""pslq_gate.py — (docstring-only edit) shared PSLQ closure harness')
    assert doc != src
    eng2 = tmp_path / 'pslq_gate_doc.py'
    eng2.write_text(doc)
    rc, text = run('pslq_calpha.py', base + ['--engine', str(eng2), '--out', str(tmp_path / 'c6b.json')], tmp_path)
    assert rc == 3 and 'not in KNOWN_ENGINE_FILES; protocol units match the pin' in text and 'HIT at' not in text
    rc, text = run('pslq_calpha.py', base + ['--engine', str(eng2), '--allow-unlisted-engine', '--out', str(tmp_path / 'c6c.json')], tmp_path)
    assert rc == 0 and 'WARNING: UNLISTED ENGINE FILE' in text and '(UNLISTED)' in text and 'PIN OK' in text
    j = json.load(open(tmp_path / 'c6c.json'))
    assert j['producer']['engine']['label'] == 'UNLISTED' and j['producer']['engine']['unlisted_allowed'] is True
    rc, text = run('pslq_calpha.py', base + ['--engine', str(eng2), '--strict-engine', '--out', str(tmp_path / 'c6d.json')], tmp_path)
    assert rc == 2 and 'unrecognized arguments: --strict-engine' in text


def test_control_monkeypatched_engine_after_the_units_is_refused_by_default(tmp_path):
    """The bypass the verifier found: a module-level statement AFTER the 15 protocol units (globals()[...] = fake)
    is invisible to the token pin; the listed-file gate refuses it by default.  Under --allow-unlisted-engine the
    run proceeds and the fake (a NULL turned into (1, 1, ..., 1) at rung 1e4, reverify stubbed ok) closes a target
    that is a true NULL on the listed engine -- that is the stated exception, and the warning line names it.  The
    target is the CY3 c_{5/4} string under --family CY4 (a true NULL in the CY4 ring), so the family's recorded
    closed form is switched off (--closed none)."""
    src = open(os.path.join(PKG, '..', 'lockpick', 'pslq_gate.py')).read()
    patched = src + '''

_orig_tps = two_prec_stable
def _tps(target_str, names, members, dps_pair=DEFAULT['dps_pair'], maxcoeff=DEFAULT['maxcoeff'], **kw):
    v = _orig_tps(target_str, names, members, dps_pair=dps_pair, maxcoeff=maxcoeff, **kw)
    if v is None and int(maxcoeff) == 10**4:
        return tuple([1] + [1] * len(names))
    return v
def _rev(hit, members, **kw):
    return {'relresid': '0.0', 'ok': True}
globals()['two_prec_stable'] = _tps
globals()['reverify'] = _rev
'''
    eng = tmp_path / 'pslq_gate_monkeypatch.py'
    eng.write_text(patched)
    assert P.protocol_sha(str(eng)) == P.PSLQ_GATE_PROTOCOL_SHA256          # the token pin alone cannot see it
    base = ['--family', 'CY4', '--c', '3/2=' + RECORD['c54'], '--form', 'log', '--dps-pair', '120,200', '--closed', 'none']
    rc, text = run('pslq_calpha.py', base + ['--engine', str(eng), '--out', str(tmp_path / 'mp.json')], tmp_path)
    assert rc == 3 and 'not in KNOWN_ENGINE_FILES; protocol units match the pin' in text and 'HIT at' not in text
    rc, text = run('pslq_calpha.py', base + ['--engine', str(eng), '--allow-unlisted-engine', '--out', str(tmp_path / 'mp2.json')], tmp_path)
    assert rc == 0 and 'WARNING: UNLISTED ENGINE FILE' in text and 'HIT at rung 10000' in text   # the stated exception, loud
    rc, text = run('pslq_calpha.py', base + ['--out', str(tmp_path / 'null.json')], tmp_path)
    assert rc == 1 and 'NULL-with-floor' in text and 'HIT at' not in text          # the listed engine: a true NULL


def test_control_comment_byte_tamper_is_refused_by_default(tmp_path):
    src = open(os.path.join(PKG, '..', 'lockpick', 'pslq_gate.py')).read()
    bad = src.replace('            if r[0] == 0:   # internal relation among members shadows the target\n',
                      '            if r[0] == 0:   # internal relation among members shadows the target (comment-byte tamper)\n')
    assert bad != src
    eng = tmp_path / 'pslq_gate_comment.py'
    eng.write_text(bad)
    assert P.protocol_sha(str(eng)) == P.PSLQ_GATE_PROTOCOL_SHA256          # comments are stripped from the token pin
    rc, text = run('pslq_calpha.py', ['--family', 'CY3', '--from-routeB', fx('ROUTEB_CY3_hi.json'), '--dps-pair', '120,200',
                                      '--form', 'log', '--engine', str(eng), '--out', str(tmp_path / 'cb.json')], tmp_path)
    assert rc == 3 and 'REFUSED: engine file sha256' in text and 'protocol units match the pin' in text and 'HIT at' not in text


def test_control_unparseable_engine_is_refused_by_name(tmp_path):
    eng = tmp_path / 'broken.py'
    eng.write_text('def x(:\n')
    rc, text = run('pslq_calpha.py', ['--family', 'CY3', '--from-routeB', fx('ROUTEB_CY3_hi.json'), '--form', 'log',
                                      '--engine', str(eng), '--out', str(tmp_path / 'br.json')], tmp_path)
    assert rc == 3 and 'REFUSED: engine' in text and 'unparseable engine file (SyntaxError' in text and 'Traceback' not in text


def cy4_planted_string(digits=320):
    """-sqrt(5)/(100 pi) computed at run time (never pasted): a value in the CY4 predicted ring."""
    import mpmath as mp
    with mp.workdps(digits + 20):
        return mp.nstr(-mp.sqrt(5) / (100 * mp.pi), digits, strip_zeros=False)


def test_three_member_ring_controls_pass_and_close(tmp_path):
    """A 3-member --ring basket: the planted control's third index used to coincide with the second (every 3-member
    basket failed its own control, rc 2, no verdict); it now plants three distinct members.  The targets are synthetic
    values in the CY4 ring, so the family's recorded closed form is switched off (--closed none)."""
    t = cy4_planted_string()
    r3 = tmp_path / 'ring3_log.json'
    json.dump({"log": ["log_pi", "log2", "log5"]}, open(r3, 'w'))          # log|c| = -log_pi - 2 log2 - (3/2) log5
    rc, text = run('pslq_calpha.py', ['--family', 'CY4', '--c', '3/2=' + t, '--ring', str(r3), '--form', 'log', '--closed', 'none',
                                      '--out', str(tmp_path / 'r3.json')], tmp_path)
    assert rc == 0, text
    j = json.load(open(tmp_path / 'r3.json'))
    s = j['scans'][0]
    assert s['names'] == ['log_pi', 'log2', 'log5'] and tuple(s['HIT']['vector']) == (2, 2, 4, 3) and s['verdict'] == 'CLOSED'
    assert all(c['ok'] and all(d['ok'] for d in c['detail']) for c in j['controls']['log'].values())
    assert 'control planted (3*log2 - 7*log5 + 2*log_pi)/5: PASS' in text
    r3a = tmp_path / 'ring3_add.json'
    json.dump({"additive": {"pi_power": "1", "members": ["1", "sqrt5", "sqrt3"]}}, open(r3a, 'w'))   # c*pi = -sqrt5/100
    rc, text = run('pslq_calpha.py', ['--family', 'CY4', '--c', '3/2=' + t, '--ring', str(r3a), '--form', 'additive', '--closed', 'none',
                                      '--out', str(tmp_path / 'r3a.json')], tmp_path)
    assert rc == 0, text
    s = json.load(open(tmp_path / 'r3a.json'))['scans'][0]
    assert tuple(s['HIT']['vector']) == (100, 0, 1, 0) and 'control planted (3*sqrt5 - 7*sqrt3 + 2*1)/5: PASS' in text
    r3n = tmp_path / 'ring3_null.json'
    json.dump({"log": ["1", "log_pi", "log2"]}, open(r3n, 'w'))              # log5 absent: a true NULL, controls green
    rc, text = run('pslq_calpha.py', ['--family', 'CY4', '--c', '3/2=' + t, '--ring', str(r3n), '--form', 'log', '--closed', 'none',
                                      '--out', str(tmp_path / 'r3n.json')], tmp_path)
    assert rc == 1 and 'NULL-with-floor' in text and 'CONTROLS FAILED' not in text
    assert 'control planted (3*log_pi - 7*log2 + 2*1)/5: PASS' in text


def test_control_wrong_closed_form_is_a_named_mismatch(tmp_path):
    wrong = tmp_path / 'wrong.json'
    json.dump({"5/4": {"sign": -1, "rational": "1/2", "pi": "-1/2", "radicals": {"2": "-1/2"}, "gamma": {"1/4": -2},
                       "display": "WRONG"}}, open(wrong, 'w'))
    rc, text = run('pslq_calpha.py', ['--family', 'CY3', '--from-routeB', fx('ROUTEB_CY3_hi.json'), '--dps-pair', '120,200',
                                      '--form', 'log', '--closed', str(wrong), '--out', str(tmp_path / 'c2.json')], tmp_path)
    assert rc == 1 and 'EXPECTATION MISMATCH: expected (2, 0, 1, 3, 0, 0, 0, 4)' in text and 'got (2, 0, 1, 1, 0, 0, 0, 4)' in text


def test_control_failed_controls_give_no_verdict(tmp_path):
    rc, text = run('pslq_calpha.py', ['--family', 'CY3', '--from-routeB', fx('ROUTEB_CY3_hi.json'), '--dps-pair', '120,200',
                                      '--form', 'log', '--maxcoeff', '5', '--ladder', '5', '--out', str(tmp_path / 'c7.json')], tmp_path)
    assert rc == 2 and 'CONTROLS FAILED' in text and 'no verdict' in text and 'HIT at' not in text

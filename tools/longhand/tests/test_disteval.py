# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""The battery of Longhand's disteval route (tools/longhand/disteval; the worked example's records under tools/longhand/examples/ndpent_top).
The worked-example leg: the acceptance check (--stage check, older spelling --stage gate) over the shipped records must reproduce the record (the ndpent top-sector scalar at the
Euclidean point s12=-7 s23=-9 s34=-13 s45=-10 s15=-5 mm3=-12 mm4=-8; 332 sectors, 968 kernels) -- the conservative digit class 7, the
per-order table 9/7/7/7/7 (eps^-4 .. eps^0), 0 B-C lines over 3 sigma, the A-consistency MARGINAL verdict at max sigma 3.0519681777552745,
the PLANTED-FAIL control raised (sigma vs A 9.99707344790171), the FOREIGN controls within 3x the summed bars (max sigmas 1.717724232013687 and
2.9485284820009676) -- every figure read from GATE_record.json AND asserted as the record's literal beside it; the compare rule on a synthetic pair;
the chunk split on the package JSON (332 sectors / 968 kernels / the modal histogram {3: 122, 2: 112, 4: 72, 5: 16, 1: 10}; 28 chunks of 12 sectors);
the per-chunk epsabs rule ('1.706e-10' == the record's stage-B setting); the memory basis (W 8 PROBE 377241600 B == the fence of record; W 12/16/24
per the record); the refusals by name (a package pin mismatch rc 3, a missing fixture rc 4, a bad --point rc 2, --unpinned-package's loud line);
the dry (stand-in) stage B with its checkpoints and --resume; the controls (a planted digit in RESULT_B re-pinned -> the gate FAILS by name; a
smoke shifted outside 3x the bars -> foreign_ok False).  Fixtures sha256-pinned in PINS.json; the compiled package is NOT a fixture (the
pins in PACKAGE_PINS.json; the on-host check is skipped by name when the package dir is absent).
The family legs: --family ndpent_top explicit == the default on every tier (masked captures identical); a SYNTHETIC family (the fixture
package JSONs renamed to ndpent_m3_n123) through the chunk split, the dry stage B with checkpoints and --resume, the sums-key path of the
assembly, the gate's smoke parse (a smoke of another family FAILS the FOREIGN control by name); --pins FILE (match OK, one byte planted rc 3,
the pins file missing rc 4; the default pins unchanged); --emit-pins; the --name alias on the helper CLIs; a --family naming no files rc 2.
The refusal legs (one per refusal, each asserting the named line AND the rc, never a traceback): a --family of the wrong form (empty, a
separator, a dot, a path) rc 2 on every arm and every helper CLI before anything is read or written; a sum file whose 'integrals' is not
[NAME_integral], an integral file whose 'name' is not NAME_integral, a sum file that is not JSON: rc 2 on the stage line, --check-package,
--emit-pins (nothing written) and stage_chunks.py; a --pins file of the wrong shape (no "files", not JSON, a list, "files" not an object, a
value not 64-hex, an absolute key, "family" not a string): rc 2 on --check-package and the stage line, --unpinned-package no bypass; a
--generate-receipt file missing: rc 4, nothing written; the helper CLIs without the family of the package (stage_chunks.py, assemble) and a
chunks dir without its manifest: rc 2 by name.
Every subprocess carries the interpreter's user site (PYTHONUSERBASE) so a scratch HOME never hides the user-site pySecDec."""
import hashlib, json, math, os, shutil, site, subprocess, sys, types

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                       # tools/longhand, the package root
PKG = os.path.join(ROOT, 'disteval')               # the disteval member: the CLI and its helper modules
FIX = os.path.join(ROOT, 'examples', 'ndpent_top')  # the worked example's records, sha256-pinned in PINS.json
EX_REL = os.path.join('examples', 'ndpent_top')    # the same, relative to a copied package root
CLI_REL = os.path.join('disteval', 'pysecdec_point.py')
PY = sys.executable
CLI = os.path.join(PKG, 'pysecdec_point.py')
sys.path.insert(0, PKG)
import stage_chunks, assemble_compare, gate as gate_mod, memory_basis  # noqa: E402

RECORD = {'digit_class': 7, 'per_order': [9, 7, 7, 7, 7], 'A_max_sigma': 3.0519681777552745, 'planted_sigma_A': 9.99707344790171, 'foreign_max': [1.717724232013687, 2.9485284820009676],
          'eps0_B': 0.8246341248952935, 'eps0_B_err': 7.344564291155832e-08, 'n_sectors': 332, 'n_kernels': 968, 'hist': {3: 122, 2: 112, 4: 72, 5: 16, 1: 10}, 'typical_sector': 151,
          'n_chunks_B': 28, 'epsabs_B': '1.706e-10', 'memmax_W8_probe': 377241600, 'memmax': {8: 377241600, 12: 222425088, 16: 256225280, 24: 323829760}, 'worker_B': 6500352, 'driver_B': 93089792}
PACKAGE_PINS_SHA16 = 'ce7513935285834d'   # the fixture family's PACKAGE_PINS.json, unchanged by the --family / --pins succession
SYN = 'ndpent_m3_n123'                     # the synthetic family name of the family legs


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def fx(n):
    return os.path.join(FIX, n)


def base_env():
    """the subprocess env: PATH, HOME when set, and the user site (PYTHONUSERBASE when set, else the interpreter's own user base read at run
    time) -- a scratch HOME alone would hide the user-site pySecDec from the stand-in check; PYTHONDONTWRITEBYTECODE=1."""
    env = {'PATH': os.environ.get('PATH', '/usr/bin:/bin'), 'PYTHONDONTWRITEBYTECODE': '1'}
    if 'HOME' in os.environ:
        env['HOME'] = os.environ['HOME']
    ub = os.environ.get('PYTHONUSERBASE') or site.getuserbase()
    if ub:
        env['PYTHONUSERBASE'] = ub
    return env


def run(args, cwd=None, env_extra=None):
    env = base_env()
    env.update(env_extra or {})
    p = subprocess.run([PY, CLI] + args, cwd=cwd or HERE, capture_output=True, text=True, env=env)
    return p.returncode, p.stdout + p.stderr


def run_py(script, args, cwd=None):
    p = subprocess.run([PY, script] + args, cwd=cwd or HERE, capture_output=True, text=True, env=base_env())
    return p.returncode, p.stdout + p.stderr


def mask(txt):
    """the run-dependent tokens of a capture: date -u stamps, sha16 tokens of freshly written files, GNU-time walls and rss figures."""
    import re
    txt = re.sub(r'\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ', '<STAMP>', txt)
    txt = re.sub(r'\b[0-9a-f]{16}\b', '<SHA16>', txt)
    txt = re.sub(r'wall=\S+ maxrss_kb=\S+', 'wall=<W> maxrss_kb=<R>', txt)
    return txt


def have_pysecdec():
    try:
        import pySecDec.disteval, pySecDecContrib  # noqa: F401
        return os.access(os.path.join(os.path.dirname(pySecDecContrib.__file__), 'bin', 'pysecdec_cpuworker'), os.X_OK)
    except Exception:
        return False


def fake_package(root, with_so=True):
    """a package dir for the split / dry legs: the two package JSONs of record + zero-length .so stand-ins + an empty coefficients dir."""
    os.makedirs(os.path.join(root, 'coefficients'), exist_ok=True)
    shutil.copyfile(fx('package_ndpent_top.json'), os.path.join(root, 'ndpent_top.json'))
    shutil.copyfile(fx('package_ndpent_top_integral.json'), os.path.join(root, 'ndpent_top_integral.json'))
    open(os.path.join(root, 'coefficients', 'ndpent_top_integral_coefficient0.txt'), 'w').write('1\n')
    if with_so:
        for so in ('ndpent_top_integral.so', 'builtin.so'):
            open(os.path.join(root, so), 'wb').write(b'not-a-library\n')
    return root


def synthetic_package(root, name=SYN):
    """the SYNTHETIC family: the fixture package JSONs with every name re-cut to NAME (the sum file's name / integrals / sums key / coefficient,
    the integral's name); the kernels, orders and the prefactor untouched; .so stand-ins + the coefficient file as fake_package."""
    os.makedirs(os.path.join(root, 'coefficients'), exist_ok=True)
    S = json.load(open(fx('package_ndpent_top.json'))); I = json.load(open(fx('package_ndpent_top_integral.json')))
    S['name'] = name; S['integrals'] = [name + '_integral']
    S['sums'] = {name: [dict(e, integral=name + '_integral', coefficient=name + '_integral_coefficient0.txt') for e in S['sums']['ndpent_top']]}
    I['name'] = name + '_integral'
    json.dump(S, open(os.path.join(root, name + '.json'), 'w'), indent=2); json.dump(I, open(os.path.join(root, name + '_integral.json'), 'w'), indent=2)
    open(os.path.join(root, 'coefficients', name + '_integral_coefficient0.txt'), 'w').write('1\n')
    for so in (name + '_integral.so', 'builtin.so'):
        open(os.path.join(root, so), 'wb').write(b'not-a-library\n')
    return root


def synthetic_smoke(out, name=SYN):
    """the fixture smoke with its sums / integrals keys re-cut to the synthetic family (the rows untouched): the dry stand-in result + a FOREIGN smoke."""
    sm = json.load(open(fx('smoke_build_host.json')))
    sm['sums'] = {name: sm['sums']['ndpent_top']}; sm['integrals'] = {name + '_integral': sm['integrals']['ndpent_top_integral']}
    json.dump(sm, open(out, 'w'), indent=1)
    return out


# ---------------------------------------------------------------- pins
def test_fixture_pins():
    pins = json.load(open(fx('PINS.json')))
    files = sorted(f for f in os.listdir(FIX) if f != 'PINS.json')
    assert sorted(pins) == files
    for f in files:
        assert sha(fx(f)) == pins[f], f
    assert len(pins) == 17


def test_package_pins_name_the_fixture_family():
    pk = json.load(open(os.path.join(PKG, 'PACKAGE_PINS.json')))
    assert pk['family'] == 'ndpent_top' and pk['integral'] == 'ndpent_top_integral'
    assert set(pk['files']) == {'builtin.so', 'ndpent_top.json', 'ndpent_top_integral.json', 'ndpent_top_integral.so', 'coefficients/ndpent_top_integral_coefficient0.txt'}
    assert pk['files']['ndpent_top_integral.json'] == sha(fx('package_ndpent_top_integral.json'))
    assert pk['files']['ndpent_top.json'] == sha(fx('package_ndpent_top.json'))
    assert pk['n_sectors'] == RECORD['n_sectors'] and pk['n_kernels'] == RECORD['n_kernels']
    assert pk['point_of_record'] == 's12=-7 s23=-9 s34=-13 s45=-10 s15=-5 mm3=-12 mm4=-8'


def test_package_on_this_host_matches_the_pins_or_is_absent():
    d = os.environ.get('LONGHAND_DISTEVAL_PACKAGE_DIR') or os.environ.get('PYSECDEC_POINT_PACKAGE_DIR')
    if not d or not os.path.isdir(d):
        pytest.skip('no compiled package on this host (LONGHAND_DISTEVAL_PACKAGE_DIR, or the older name PYSECDEC_POINT_PACKAGE_DIR, unset): the package is not a fixture')
    rc, out = run(['--check-package', '--disteval-dir', d])
    assert rc == 0 and 'PACKAGE PINNED' in out, out


# ---------------------------------------------------------------- the gate reproduces the record
def test_gate_reproduces_the_record(tmp_path):
    out = str(tmp_path / 'GATE.json')
    rc, txt = run(['--stage', 'gate', '--fixtures', '--out', out])
    assert rc == 0, txt
    assert 'RECORD REPRODUCED: 21 checks' in txt
    g = json.load(open(out)); rec = json.load(open(fx('GATE_record.json')))
    assert g['gate'] == 'PASS' == rec['gate']
    assert g['digit_class']['min_digits_over_orders(conservative)'] == RECORD['digit_class'] == rec['digit_class']['min_digits_over_orders(conservative)']
    key = 'digits(min over the B-C pair and the B, C bars)'
    assert [r[key] for r in g['digit_class']['per_order']] == RECORD['per_order'] == [r[key] for r in rec['digit_class']['per_order']]
    assert [r['regulator_powers'] for r in g['digit_class']['per_order']] == [[-4], [-3], [-2], [-1], [0]]
    assert g['fail_lines(B-C sigma > 3 on an order)'] == [] and g['A_consistency']['max_sigma_A_vs_B_or_C'] == RECORD['A_max_sigma'] == rec['A_consistency']['max_sigma_A_vs_B_or_C']
    assert g['A_consistency']['verdict'].startswith('MARGINAL (2 line(s) over 3 sigma')
    assert g['controls']['planted_fail']['raised_fail_line'] is True and g['controls']['planted_fail']['sigma_vs_A'] == RECORD['planted_sigma_A']
    assert g['controls']['foreign_ok'] is True and [f['max_sigma'] for f in g['controls']['foreign']] == RECORD['foreign_max']
    for r, rr in zip(g['digit_class']['per_order'], rec['digit_class']['per_order']):
        for k in ('A-B_sigma', 'A-C_sigma', 'B-C_sigma', 'A_err', 'B_err', 'C_err', 'A_value', 'B_value', 'C_value'):
            assert r[k] == rr[k], k
    e0 = [r for r in g['digit_class']['per_order'] if r['regulator_powers'] == [0]][0]
    assert e0['B_value'][0] == RECORD['eps0_B'] and e0['B_err'] == RECORD['eps0_B_err']
    assert g['package_pinned'] is True and g['leaf_at_end']['oom_kill'] == '0' and g['leaf_at_end']['memory.max'] == '377241600'


def test_stage_check_and_example_are_the_acceptance_check_aliases(tmp_path):
    """--stage check --example is the acceptance check over the worked example: the same computation as the older spelling --stage gate --fixtures
    (masked captures identical, the same GATE json but for the stamp), and the refusal line names both spellings."""
    o1, o2 = str(tmp_path / 'G1.json'), str(tmp_path / 'G2.json')
    rc1, t1 = run(['--stage', 'gate', '--fixtures', '--out', o1]); rc2, t2 = run(['--stage', 'check', '--example', '--out', o2])
    assert rc1 == rc2 == 0 and 'RECORD REPRODUCED: 21 checks' in t2, t2
    assert mask(t1).replace(o1, '<OUT>') == mask(t2).replace(o2, '<OUT>')
    g1, g2 = json.load(open(o1)), json.load(open(o2))
    g1.pop('stamp_utc'); g2.pop('stamp_utc')
    assert g1 == g2 and g2['gate'] == 'PASS'
    rc, txt = run(['--stage', 'check']); assert rc == 2 and '--example/--fixtures or --work' in txt and 'Traceback' not in txt, txt
    rc, txt = run(['--stage', 'check', '--fixtures', '--out', str(tmp_path / 'G3.json')]); assert rc == 0      # mixed spellings


def test_gate_direct_over_the_fixture_results_matches_the_cli(tmp_path):
    out = str(tmp_path / 'GATE_direct.json')
    res = {s: fx('RESULT_%s.json' % s) for s in 'ABC'}
    g = gate_mod.run_gate(res, {s: 'rc=0' for s in 'ABC'}, [('build host', fx('smoke_build_host.json')), ('second host', fx('smoke_second_host.json'))], gate_mod.read_stamp(fx('LEAF_STAMP_record.txt')), out)
    assert g['gate'] == 'PASS' and g['digit_class']['min_digits_over_orders(conservative)'] == 7


# ---------------------------------------------------------------- the controls
def test_planted_fail_control_line():
    rc, txt = run(['--planted-fail-control'])
    assert rc == 0 and 'FAIL line RAISED' in txt and '9.99707344790171' in txt, txt


def test_control_planted_digit_in_result_b_repinned_fails_the_gate_by_name(tmp_path):
    """a copy of the tool with RESULT_B's eps^0 value shifted by 1e-3 relative and PINS re-pinned: the gate FAILS by name (rc 1) on the B-C line at eps^0."""
    T = str(tmp_path / 'tool'); shutil.copytree(ROOT, T, ignore=shutil.ignore_patterns('__pycache__', '.pytest_cache'))
    fb = os.path.join(T, EX_REL, 'RESULT_B.json'); j = json.load(open(fb))
    for o in j['orders']:
        if o['regulator_powers'] == [0]:
            o['value_re'] = o['value_re'] * (1 + 1e-3)
    json.dump(j, open(fb, 'w'), indent=1)
    pins = json.load(open(os.path.join(T, EX_REL, 'PINS.json'))); pins['RESULT_B.json'] = sha(fb); json.dump(pins, open(os.path.join(T, EX_REL, 'PINS.json'), 'w'), indent=1)
    p = subprocess.run([PY, os.path.join(T, CLI_REL), '--stage', 'gate', '--fixtures', '--out', str(tmp_path / 'G.json')], capture_output=True, text=True, env=base_env())
    assert p.returncode == 1, p.stdout + p.stderr
    g = json.load(open(str(tmp_path / 'G.json')))
    assert g['gate'] == 'FAIL' and g['fail_lines(B-C sigma > 3 on an order)'][0]['regulator_powers'] == [0] and g['fail_lines(B-C sigma > 3 on an order)'][0]['sigma'] > 3
    assert 'FAIL line: pair B-C order [0]' in p.stdout and 'RECORD NOT REPRODUCED' in p.stdout


def test_control_planted_digit_without_repin_is_a_pin_mismatch(tmp_path):
    T = str(tmp_path / 'tool'); shutil.copytree(ROOT, T, ignore=shutil.ignore_patterns('__pycache__', '.pytest_cache'))
    fb = os.path.join(T, EX_REL, 'RESULT_B.json'); j = json.load(open(fb)); j['orders'][0]['value_re'] *= (1 + 1e-9); json.dump(j, open(fb, 'w'), indent=1)
    p = subprocess.run([PY, os.path.join(T, CLI_REL), '--stage', 'gate', '--fixtures', '--out', str(tmp_path / 'G.json')], capture_output=True, text=True, env=base_env())
    assert p.returncode == 3 and 'fixture pin mismatch: RESULT_B.json' in p.stdout, p.stdout + p.stderr


def test_foreign_control_rows_and_a_shifted_smoke_fails(tmp_path):
    rec = json.load(open(fx('GATE_record.json')))
    O = {s: gate_mod.orders_of(json.load(open(fx('RESULT_%s.json' % s)))) for s in 'ABC'}
    S = gate_mod.smoke_orders(fx('smoke_build_host.json'))
    pr = gate_mod.pair(O['A'], S)
    assert all(r['within_3x_sum_of_bars'] for r in pr.values()) and max(r['sigma'] for r in pr.values()) == RECORD['foreign_max'][0]
    assert {str(list(k)): v for k, v in pr.items()} == rec['controls']['foreign'][0]['rows']
    # a smoke shifted by 10x the summed bars on eps^0 -> outside 3x -> foreign_ok False -> the gate FAILS by name
    sm = json.load(open(fx('smoke_build_host.json')))
    for row in sm['sums']['ndpent_top']:
        if row[0] == [0]:
            row[1][0] += 10 * (row[2][0] + O['A'][(0,)][1])
    bad = str(tmp_path / 'smoke_shifted.json'); json.dump(sm, open(bad, 'w'))
    g = gate_mod.run_gate({s: fx('RESULT_%s.json' % s) for s in 'ABC'}, {s: 'rc=0' for s in 'ABC'}, [('shifted', bad)], {}, str(tmp_path / 'G.json'))
    assert g['gate'] == 'FAIL' and g['controls']['foreign_ok'] is False and g['controls']['foreign'][0]['rows']['[0]']['within_3x_sum_of_bars'] is False


def test_control_no_foreign_smoke_is_never_a_pass(tmp_path):
    g = gate_mod.run_gate({s: fx('RESULT_%s.json' % s) for s in 'ABC'}, {s: 'rc=0' for s in 'ABC'}, [], {}, str(tmp_path / 'G.json'))
    assert g['gate'] == 'FAIL' and g['controls']['foreign_ok'] is False and g['controls']['foreign'] == []


# ---------------------------------------------------------------- the compare rule
def test_compare_rule_on_a_synthetic_pair(tmp_path):
    def result(name, v, e):
        return {'complete': True, 'orders': [{'regulator_powers': [0], 'value_re': v, 'value_im': 0.0, 'err_re': e, 'err_im': 0.0, 'bar_digits': math.floor(-math.log10(e / abs(v)))}]}
    a, b = str(tmp_path / 'A.json'), str(tmp_path / 'B.json')
    json.dump(result('A', 1.0, 1e-3), open(a, 'w')); json.dump(result('B', 1.0015, 1e-3), open(b, 'w'))
    out = str(tmp_path / 'C.json')
    assemble_compare.compare(types.SimpleNamespace(results=['A=' + a, 'B=' + b], out=out))
    c = json.load(open(out)); r = c['pairs'][0]['rows'][0]
    assert r['sigma'] == pytest.approx(0.0015 / math.sqrt(2e-6)) and r['agreed_digits'] == math.floor(-math.log10(0.0015 / 1.0015)) == 2
    assert r['within_sum_of_bars'] is True and r['bar_digits_min'] == 3 and c['fail_lines'] == []   # bar digits floor(-log10(1e-3 / 1.0)) = 3; the conservative count = min(agreed 2, bars 3) = 2
    # the boundary: a difference of exactly the summed bars in decimal (1.002 - 1.0 vs 0.002) lands on the wrong side of the double comparison; the rule is the plain <=
    json.dump(result('B', 1.002, 1e-3), open(b, 'w')); assemble_compare.compare(types.SimpleNamespace(results=['A=' + a, 'B=' + b], out=out))
    r = json.load(open(out))['pairs'][0]['rows'][0]; assert r['within_sum_of_bars'] == (abs(1.002 - 1.0) <= 2e-3) and r['sigma'] == pytest.approx(0.002 / math.sqrt(2e-6))
    assert c['min_digits_over_all_orders_and_pairs(the conservative digit count)'] == 2
    json.dump(result('B', 1.01, 1e-3), open(b, 'w')); assemble_compare.compare(types.SimpleNamespace(results=['A=' + a, 'B=' + b], out=out))
    c = json.load(open(out)); assert len(c['fail_lines']) == 1 and c['fail_lines'][0]['sigma'] > 3


def test_compare_over_the_fixture_results_reproduces_the_record(tmp_path):
    out = str(tmp_path / 'COMPARE.json')
    rc, txt = run(['--stage', 'compare', '--results'] + ['%s=%s' % (s, fx('RESULT_%s.json' % s)) for s in 'ABC'] + ['--out', out])
    assert rc == 0, txt
    c = json.load(open(out)); rec = json.load(open(fx('COMPARE_record.json')))
    assert c['min_digits_over_all_orders_and_pairs(the conservative digit count)'] == rec['min_digits_over_all_orders_and_pairs(the conservative digit count)'] == 3
    assert c['fail_lines'] == rec['fail_lines'] and len(c['fail_lines']) == 2 and {f['pair'][0] for f in c['fail_lines']} == {'A'}
    for pc, pr in zip(c['pairs'], rec['pairs']):
        assert pc['pair'] == pr['pair']
        for rc_, rr in zip(pc['rows'], pr['rows']):
            assert rc_['sigma'] == rr['sigma'] and rc_['agreed_digits'] == rr['agreed_digits'] and rc_['abs_diff'] == rr['abs_diff']


# ---------------------------------------------------------------- the chunk split
def test_chunk_split_on_the_package_json(tmp_path):
    D = fake_package(str(tmp_path / 'pkg'))
    S, I = stage_chunks.load(D)
    secs = stage_chunks.sectors_in_order(I)
    assert len(secs) == RECORD['n_sectors'] and len(I['kernels']) == RECORD['n_kernels']
    s, mode, ncand, hist = stage_chunks.choose_typical(I)
    assert (s, mode, ncand) == (RECORD['typical_sector'], 3, 122) and hist == RECORD['hist']
    p = subprocess.run([PY, os.path.join(PKG, 'stage_chunks.py'), '--disteval-dir', D, '--out', str(tmp_path / 'chunks'), '--chunk-sectors', '12', '--stage', 'B'], capture_output=True, text=True, env=base_env())
    assert p.returncode == 0, p.stdout + p.stderr
    man = json.load(open(str(tmp_path / 'chunks' / 'CHUNKS_B.json')))
    assert man['n_chunks'] == RECORD['n_chunks_B'] and man['kernels_cover_package'] is True and man['n_kernels_in_chunks'] == RECORD['n_kernels']
    rec = json.load(open(fx('CHUNKS_B_record.json')))
    assert [c['sectors'] for c in man['chunks']] == [c['sectors'] for c in rec['chunks']]
    assert [c['n_kernels'] for c in man['chunks']] == [c['n_kernels'] for c in rec['chunks']]
    assert [c['integral_json_sha256'] for c in man['chunks']] == [c['integral_json_sha256'] for c in rec['chunks']]
    c0 = json.load(open(os.path.join(man['chunks'][0]['dir'], 'ndpent_top_integral.json')))
    assert all(stage_chunks.sector_of(k) in set(man['chunks'][0]['sectors']) for k in c0['kernels']) and all(stage_chunks.sector_of(k) in set(man['chunks'][0]['sectors']) for o in c0['orders'] for k in o['kernels'])
    assert os.path.islink(os.path.join(man['chunks'][0]['dir'], 'coefficients'))
    p = subprocess.run([PY, os.path.join(PKG, 'stage_chunks.py'), '--disteval-dir', D, '--choose-typical'], capture_output=True, text=True, env=base_env())
    j = json.loads(p.stdout); assert j['sector'] == 151 and j['modal_kernel_count'] == 3 and j['n_sectors'] == 332 and j['n_kernels'] == 968


def test_epsabs_per_chunk_rule_reproduces_the_record_setting():
    import pysecdec_point as cli
    assert cli.epsabs_per_chunk(fx('RESULT_A.json'), '1e-7', 28, '1e-12') == RECORD['epsabs_B'] == json.load(open(fx('SETTINGS_B_record.json')))['epsabs_per_chunk'] == json.load(open(fx('RESULT_B.json')))['settings']['epsabs_per_chunk']
    assert cli.epsabs_per_chunk(fx('RESULT_A.json'), '1e-7', 28, '1e-9') == '1.000e-09'


# ---------------------------------------------------------------- the memory basis
def test_basis_formula_reproduces_the_record():
    rec = json.load(open(fx('BASIS_record.json')))
    worker_B, driver_B, _ = memory_basis.read_basis_file(fx('BASIS_record.json'))
    assert (worker_B, driver_B) == (RECORD['worker_B'], RECORD['driver_B'])
    for W in (8, 12, 16, 24):
        e = memory_basis.basis_entry(W, worker_B, driver_B, 8)
        r = rec['per_W'][str(W)]
        for k in ('formula_memory_max_B(per_worker x W x 1.3, page_up)', 'corrected_memory_max_B(1.3 x (driver + W x per_worker), page_up)', 'memory_max_bytes', 'raise_once_bytes', 'formula_headroom_over_W_x_worker_B', 'driver_fits_in_formula_headroom', 'label', 'census_memg_gib'):
            assert e[k] == r[k], (W, k)
        assert e['memory_max_bytes'] == RECORD['memmax'][W]
        assert e['corrected_memory_max_B(1.3 x (driver + W x per_worker), page_up)'] == memory_basis.page_up(1.3 * (driver_B + W * worker_B))
    fence = json.load(open(fx('FENCE_AT_SPAWN_record.json')))
    assert int(fence['memory.max']) == RECORD['memmax_W8_probe'] == memory_basis.basis_entry(8, worker_B, driver_B, 8)['memory_max_bytes'] and fence['cpu.max'] == '800000 100000' and fence['fence_ok'] is True
    stamp = gate_mod.read_stamp(fx('LEAF_STAMP_record.txt'))
    assert int(stamp['memory.peak']) < int(stamp['memory.max']) == RECORD['memmax_W8_probe'] and stamp['oom_kill'] == '0'
    rc, txt = run(['--basis', '8'])
    assert rc == 0 and 'memory.max 377241600 B' in txt and 'PROBE' in txt and 'cpu.max = 800000 100000' in txt, txt
    rc, txt = run(['--basis', '12'])
    assert rc == 0 and 'memory.max 222425088 B' in txt and 'UNCONTENDED' in txt, txt
    assert memory_basis.page_up(4097) == 8192 and memory_basis.page_up(4096) == 4096


# ---------------------------------------------------------------- refusals by name
def test_refusal_package_pin_mismatch_rc3_and_unpinned_is_loud(tmp_path):
    D = fake_package(str(tmp_path / 'pkg'))
    rc, txt = run(['--check-package', '--disteval-dir', D])
    assert rc == 3 and 'REFUSED: package pin mismatch' in txt and 'builtin.so' in txt and 'ndpent_top_integral.so' in txt, txt
    rc, txt = run(['--check-package', '--disteval-dir', D, '--unpinned-package'])
    assert rc == 0 and 'UNPINNED PACKAGE' in txt and 'package_pinned false' in txt, txt
    rc, txt = run(['--stage', 'A', '--disteval-dir', D, '--work', str(tmp_path / 'w')])
    assert rc == 3 and 'REFUSED: package pin mismatch' in txt and not os.path.exists(str(tmp_path / 'w' / 'stage_A' / 'chunks'))


def test_refusal_missing_fixture_rc4(tmp_path):
    T = str(tmp_path / 'tool'); shutil.copytree(ROOT, T, ignore=shutil.ignore_patterns('__pycache__', '.pytest_cache', 'RESULT_C.json'))
    p = subprocess.run([PY, os.path.join(T, CLI_REL), '--stage', 'gate', '--fixtures', '--out', str(tmp_path / 'G.json')], capture_output=True, text=True, env=base_env())
    assert p.returncode == 4 and 'pinned fixture missing' in p.stdout and 'RESULT_C.json' in p.stdout, p.stdout + p.stderr


def test_refusal_bad_point_rc2(tmp_path):
    D = fake_package(str(tmp_path / 'pkg'))
    for bad in ('s12=-7 s23=-9', 's12=-7 s23=-9 s34=-13 s45=-10 s15=-5 mm3=-12 mm4=abc', 's12=-7 s23=-9 s34=-13 s45=-10 s15=-5 mm3=-12 mm4=-8 extra=1', 's12=-7 s23=-9 s34=-13 s45=-10 s15=-5 mm3=-12 mm4=-8 mm4=-8', 'nonsense'):
        rc, txt = run(['--stage', 'A', '--disteval-dir', D, '--work', str(tmp_path / 'w'), '--unpinned-package', '--point', bad, '--dry'])
        assert rc == 2 and 'REFUSED: --point' in txt, (bad, txt)


def test_refusal_usage_rc2(tmp_path):
    rc, txt = run(['--stage', 'A', '--work', str(tmp_path / 'w')]); assert rc == 2 and 'needs --disteval-dir' in txt
    rc, txt = run(['--stage', 'gate']); assert rc == 2 and '--fixtures or --work' in txt
    rc, txt = run(['--stage', 'compare']); assert rc == 2 and '--results' in txt
    rc, txt = run([]); assert rc == 2


# ---------------------------------------------------------------- the dry (stand-in) stage: checkpoints, assembly, resume
@pytest.mark.skipif(not have_pysecdec(), reason='pySecDec 1.6.6 with pySecDecContrib is not importable on this host: the stand-in check imports it')
def test_dry_stage_b_checkpoints_and_resume(tmp_path):
    D = fake_package(str(tmp_path / 'pkg'))
    W = str(tmp_path / 'work')
    rc, txt = run(['--stage', 'B', '--disteval-dir', D, '--unpinned-package', '--work', W, '--workers', '2', '--dry', '--epsabs-from', fx('RESULT_A.json')])
    assert rc == 0 and 'STAGE B DONE' in txt and 'per-chunk epsabs := 1.706e-10' in txt, txt
    R = json.load(open(os.path.join(W, 'stage_B', 'RESULT_B.json')))
    assert R['complete'] is True and R['n_chunks'] == 28 == R['n_landed'] and R['settings']['epsabs_per_chunk'] == '1.706e-10' and R['settings']['standin'] == '1'
    assert R['settings']['package']['pinned_to_the_fixture_of_record'] is False
    canned = json.load(open(fx('smoke_build_host.json')))['sums']['ndpent_top']
    e0 = [o for o in R['orders'] if o['regulator_powers'] == [0]][0]
    v = 0.0
    for _ in range(28):
        v += [r for r in canned if r[0] == [0]][0][1][0]
    assert e0['value_re'] == v and e0['err_re'] == pytest.approx(math.sqrt(28) * [r for r in canned if r[0] == [0]][0][2][0])
    for ci in range(28):
        cd = os.path.join(W, 'stage_B', 'chunks', 'chunk_%03d' % ci)
        assert os.path.exists(os.path.join(cd, 'DONE')) and open(os.path.join(cd, 'rc')).read().strip() == 'rc=0' and os.path.exists(os.path.join(cd, 'time.txt'))
    assert os.path.exists(os.path.join(W, 'stage_B', 'STAGE_DONE'))
    n_prog = len(open(os.path.join(W, 'progress.jsonl')).read().splitlines()); assert n_prog == 28
    rc, txt = run(['--stage', 'B', '--disteval-dir', D, '--unpinned-package', '--work', W, '--workers', '2', '--dry', '--resume'])
    assert rc == 0 and 'already DONE' in txt and len(open(os.path.join(W, 'progress.jsonl')).read().splitlines()) == n_prog
    # --max-chunks 1 on a fresh work dir: one chunk evaluated, the stage INCOMPLETE (rc 1 by name), then --resume finishes it
    W2 = str(tmp_path / 'work2')
    rc, txt = run(['--stage', 'B', '--disteval-dir', D, '--unpinned-package', '--work', W2, '--workers', '2', '--dry', '--max-chunks', '1'])
    assert rc == 1 and 'STOPPED by --max-chunks 1' in txt and 'INCOMPLETE' in txt
    R2 = json.load(open(os.path.join(W2, 'stage_B', 'RESULT_B.json'))); assert R2['complete'] is False and R2['n_landed'] == 1 and R2['missing'] == list(range(1, 28))
    rc, txt = run(['--stage', 'B', '--disteval-dir', D, '--unpinned-package', '--work', W2, '--workers', '2', '--dry', '--resume'])
    assert rc == 0 and 'chunk 0 DONE already' in txt and 'STAGE B DONE' in txt
    # the gate over this work dir without a smoke of the unpinned library: FAIL on foreign_ok, by name
    for s in 'AC':
        shutil.copytree(os.path.join(W2, 'stage_B'), os.path.join(W2, 'stage_%s' % s)); os.replace(os.path.join(W2, 'stage_%s' % s, 'RESULT_B.json'), os.path.join(W2, 'stage_%s' % s, 'RESULT_%s.json' % s))
    rc, txt = run(['--stage', 'gate', '--work', W2, '--out', str(tmp_path / 'G.json')])
    assert rc == 1 and 'none given and the package is unpinned' in txt and json.load(open(str(tmp_path / 'G.json')))['package_pinned'] is False


# ---------------------------------------------------------------- the family legs (--family / --pins / --emit-pins / the --name alias)
def test_family_explicit_default_equals_the_default_on_every_tier(tmp_path):
    """--family ndpent_top explicit == the default: the fixture gate, the compare, the planted control, the basis, --check-package on a fake
    package (rc 3 / the labeled unpinned line), the chunk split and the assembly: masked captures identical, the output objects equal but
    for their stamps."""
    D = fake_package(str(tmp_path / 'pkg'))
    F = ['--family', 'ndpent_top']
    tiers = [(['--stage', 'gate', '--fixtures'], 'GATE'), (['--stage', 'compare', '--results'] + ['%s=%s' % (s, fx('RESULT_%s.json' % s)) for s in 'ABC'], 'COMPARE'),
             (['--planted-fail-control'], None), (['--basis', '8'], None), (['--check-package', '--disteval-dir', D], None), (['--check-package', '--disteval-dir', D, '--unpinned-package'], None)]
    for args, obj in tiers:
        outs = []
        caps = []
        for tag, extra in (('d', []), ('f', F)):
            a = list(args) + (['--out', str(tmp_path / ('%s_%s.json' % (obj, tag)))] if obj else []) + extra
            rc, txt = run(a)
            caps.append((rc, mask(txt).replace('%s_%s.json' % (obj, tag), '%s.json' % obj) if obj else mask(txt)))
            if obj:
                outs.append(json.load(open(str(tmp_path / ('%s_%s.json' % (obj, tag))))))
        assert caps[0] == caps[1], (args, caps)
        if outs:
            for o in outs:
                o.pop('stamp_utc', None)
            assert json.dumps(outs[0], sort_keys=True) == json.dumps(outs[1], sort_keys=True), args
    # the helper CLIs: stage_chunks.py with and without --family (the manifest equal but for the stamp)
    mans = []
    for tag, extra in (('d', []), ('f', F)):
        rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', D, '--out', str(tmp_path / ('chunks_' + tag)), '--chunk-sectors', '12', '--stage', 'B'] + extra)
        assert rc == 0, txt
        m = json.load(open(str(tmp_path / ('chunks_' + tag) / 'CHUNKS_B.json'))); m.pop('stamp_utc'); mans.append(json.dumps(m, sort_keys=True).replace('chunks_' + tag, 'chunks'))
    assert mans[0] == mans[1]
    assert json.loads(mans[0])['name'] == 'ndpent_top'


def test_family_naming_no_files_is_refused_by_name_rc2(tmp_path):
    D = fake_package(str(tmp_path / 'pkg'))
    rc, txt = run(['--stage', 'A', '--disteval-dir', D, '--family', 'no_such_family', '--unpinned-package', '--work', str(tmp_path / 'w'), '--dry'])
    assert rc == 2 and 'REFUSED: --family no_such_family names no files' in txt and 'no_such_family.json missing' in txt and 'families present: ndpent_top' in txt, txt
    assert not os.path.exists(str(tmp_path / 'w' / 'stage_A'))
    rc, txt = run(['--emit-pins', str(tmp_path / 'P.json'), '--disteval-dir', D, '--family', 'no_such_family'])
    assert rc == 2 and 'names no files' in txt and not os.path.exists(str(tmp_path / 'P.json'))


def test_pins_file_per_family_match_planted_byte_and_missing(tmp_path):
    """--pins FILE: --emit-pins on the synthetic package -> --check-package --pins matches (PACKAGE PINNED, rc 0); one byte planted in the
    package -> rc 3 by name; the pins file missing -> rc 4 by name; the pins of another family (the default) -> a mismatch by name (rc 3) and
    the labeled UNPINNED run; a pinned package file missing stays rc 3 (MISSING counts as a mismatch, as before); the default PACKAGE_PINS.json unchanged."""
    S = synthetic_package(str(tmp_path / 'syn'))
    P = str(tmp_path / 'PINS_syn.json')
    rc, txt = run(['--emit-pins', P, '--disteval-dir', S, '--family', SYN])
    assert rc == 0 and 'PINS EMITTED' in txt and os.path.exists(P), txt
    pins = json.load(open(P))
    assert pins['family'] == SYN and pins['integral'] == SYN + '_integral' and pins['n_sectors'] == RECORD['n_sectors'] and pins['n_kernels'] == RECORD['n_kernels']
    assert set(pins['files']) == {'builtin.so', SYN + '.json', SYN + '_integral.json', SYN + '_integral.so', 'coefficients/' + SYN + '_integral_coefficient0.txt'}
    assert all(pins['files'][k] == sha(os.path.join(S, k)) for k in pins['files']) and pins['generate_receipt_sha256'] is None
    rc, txt = run(['--check-package', '--disteval-dir', S, '--family', SYN, '--pins', P])
    assert rc == 0 and 'PACKAGE PINNED' in txt and '5 files' in txt, txt
    # the --generate-receipt sha recorded when given
    rc, txt = run(['--emit-pins', str(tmp_path / 'PINS_syn_gr.json'), '--disteval-dir', S, '--family', SYN, '--generate-receipt', fx('GENERATE_RECEIPT_record.json')])
    assert rc == 0 and json.load(open(str(tmp_path / 'PINS_syn_gr.json')))['generate_receipt_sha256'] == sha(fx('GENERATE_RECEIPT_record.json'))
    # one byte planted in the synthetic package's sum file -> rc 3 by name
    sp = os.path.join(S, SYN + '.json'); open(sp, 'a').write(' ')
    rc, txt = run(['--check-package', '--disteval-dir', S, '--family', SYN, '--pins', P])
    assert rc == 3 and 'REFUSED: package pin mismatch' in txt and SYN + '.json' in txt and 'PINS_syn.json' in txt, txt
    rc, txt = run(['--check-package', '--disteval-dir', S, '--family', SYN, '--pins', P, '--unpinned-package'])
    assert rc == 0 and 'UNPINNED PACKAGE' in txt and 'package_pinned false' in txt
    # the pins file missing -> rc 4 by name
    rc, txt = run(['--check-package', '--disteval-dir', S, '--family', SYN, '--pins', str(tmp_path / 'PINS_absent.json')])
    assert rc == 4 and 'REFUSED: pins file missing' in txt and 'PINS_absent.json' in txt, txt
    # the default pins (the fixture family) against the synthetic package -> a mismatch naming the family (rc 3); --unpinned-package runs labeled
    rc, txt = run(['--check-package', '--disteval-dir', S, '--family', SYN])
    assert rc == 3 and 'pins family ndpent_top != --family %s' % SYN in txt, txt
    rc, txt = run(['--check-package', '--disteval-dir', S, '--family', SYN, '--unpinned-package'])
    assert rc == 0 and 'UNPINNED PACKAGE' in txt and 'pins family ndpent_top != --family %s' % SYN in txt
    # a pinned package file missing stays rc 3 (as before: MISSING is a mismatch)
    os.replace(os.path.join(S, 'builtin.so'), os.path.join(S, 'builtin.so.aside'))
    rc, txt = run(['--check-package', '--disteval-dir', S, '--family', SYN, '--pins', P])
    assert rc == 3 and 'MISSING' in txt and 'builtin.so' in txt
    # the default pins file unchanged by the succession
    assert sha(os.path.join(PKG, 'PACKAGE_PINS.json'))[:16] == PACKAGE_PINS_SHA16


def test_name_alias_on_the_helper_clis(tmp_path):
    """--name (the reference form's spelling) == --family on stage_chunks.py, assemble_compare.py assemble and gate.py."""
    S = synthetic_package(str(tmp_path / 'syn'))
    for flag in ('--family', '--name'):
        rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', S, '--out', str(tmp_path / ('chunks' + flag)), '--chunk-sectors', '12', '--stage', 'B', flag, SYN])
        assert rc == 0 and 'CHUNKS B: 28 chunks, 332 sectors, 968 kernels' in txt, txt
        m = json.load(open(str(tmp_path / ('chunks' + flag) / 'CHUNKS_B.json')))
        assert m['name'] == SYN and set(m['source_sha256']) == {SYN + '.json', SYN + '_integral.json', SYN + '_integral.so', 'builtin.so'} and m['kernels_cover_package'] is True
        c0 = m['chunks'][0]['dir']
        assert os.path.exists(os.path.join(c0, SYN + '.json')) and os.path.exists(os.path.join(c0, SYN + '_integral.json')) and os.path.exists(os.path.join(c0, SYN + '_integral.so'))
        rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', S, '--choose-typical', flag, SYN])
        assert rc == 0 and json.loads(txt)['sector'] == RECORD['typical_sector']
    # the chunk split without the name on the synthetic package: the fixture family's files are absent -> REFUSED by name, rc 2, no traceback
    rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', S, '--out', str(tmp_path / 'chunks_none'), '--chunk-sectors', '12', '--stage', 'B'])
    assert rc == 2 and 'REFUSED: --family ndpent_top' in txt and 'ndpent_top.json missing' in txt and 'families present: %s' % SYN in txt and 'Traceback' not in txt, txt
    assert not os.path.exists(str(tmp_path / 'chunks_none'))
    sm = synthetic_smoke(str(tmp_path / 'smoke_syn.json'))
    # assemble_compare.py assemble --name over a hand-made chunk dir with DONE markers: the sums key read by the name
    cd = str(tmp_path / 'chunks--family')
    man = json.load(open(os.path.join(cd, 'CHUNKS_B.json')))
    for c in man['chunks'][:2]:
        shutil.copyfile(sm, os.path.join(c['dir'], 'result.json')); open(os.path.join(c['dir'], 'DONE'), 'w').write('x\n')
    for flag in ('--family', '--name'):
        rc, txt = run_py(os.path.join(PKG, 'assemble_compare.py'), ['assemble', '--chunks-dir', cd, '--stage', 'B', '--out', str(tmp_path / ('R' + flag + '.json')), flag, SYN])
        assert rc == 3 and 'landed 2/28' in txt, txt
        R = json.load(open(str(tmp_path / ('R' + flag + '.json'))))
        assert R['name'] == SYN and R['n_landed'] == 2 and [o for o in R['orders'] if o['regulator_powers'] == [0]][0]['value_re'] == 2 * json.load(open(sm))['sums'][SYN][4][1][0]
    rc, txt = run_py(os.path.join(PKG, 'assemble_compare.py'), ['assemble', '--chunks-dir', cd, '--stage', 'B', '--out', str(tmp_path / 'R_none.json')])
    # the fixture family's key is absent from the synthetic results: REFUSED by name (rc 2) naming the chunk result, the key and the keys present; no RESULT written
    assert rc == 2 and 'REFUSED: assemble --family ndpent_top' in txt and 'carries no sums.ndpent_top' in txt and 'keys: %s' % SYN in txt and 'Traceback' not in txt, txt
    assert not os.path.exists(str(tmp_path / 'R_none.json'))
    # gate.py --family / --name with the synthetic smoke over the fixture results: the smoke parsed by the family key (FOREIGN rows present)
    for flag in ('--family', '--name'):
        rc, txt = run_py(os.path.join(PKG, 'gate.py'), ['--result', 'A=' + fx('RESULT_A.json'), '--result', 'B=' + fx('RESULT_B.json'), '--result', 'C=' + fx('RESULT_C.json'), '--smoke', 'syn=' + sm, '--out', str(tmp_path / ('G' + flag + '.json')), flag, SYN])
        g = json.load(open(str(tmp_path / ('G' + flag + '.json'))))
        assert rc == 0 and g['gate'] == 'PASS' and g['controls']['foreign'][0]['max_sigma'] == RECORD['foreign_max'][0] and len(g['controls']['foreign'][0]['rows']) == 5, txt
    # the same smoke read under the fixture family: no sums.ndpent_top -> the FOREIGN control FAILS by name, never skipped
    rc, txt = run_py(os.path.join(PKG, 'gate.py'), ['--result', 'A=' + fx('RESULT_A.json'), '--result', 'B=' + fx('RESULT_B.json'), '--result', 'C=' + fx('RESULT_C.json'), '--smoke', 'syn=' + sm, '--out', str(tmp_path / 'G_wrong.json')])
    g = json.load(open(str(tmp_path / 'G_wrong.json')))
    assert rc == 1 and g['gate'] == 'FAIL' and g['controls']['foreign_ok'] is False and 'carries no sums.ndpent_top' in g['controls']['foreign'][0]['error'] and SYN in g['controls']['foreign'][0]['error'], txt


# ---------------------------------------------------------------- the refusals by name of the family / pins paths (one leg per refusal; the rc AND the line, never a traceback)
def named(rc, txt, code, *needles):
    """rc == code, every needle in the text, no traceback."""
    return rc == code and all(n in txt for n in needles) and 'Traceback' not in txt


def test_family_name_form_refused_rc2_before_anything(tmp_path):
    """a --family of the wrong form (empty; a trailing or inner separator; a path up; a dot; a space) is REFUSED by name, rc 2, on every arm of
    pysecdec_point.py and on the three helper CLIs, before any file is read or written."""
    D = fake_package(str(tmp_path / 'pkg'))
    FORM = 'is not a family name'
    for bad in ('', 'ndpent_top/', 'a/ndpent_top', '../ndpent_top', 'nd.pent', 'nd pent'):
        rc, txt = run(['--stage', 'A', '--disteval-dir', D, '--unpinned-package', '--work', str(tmp_path / 'w'), '--dry', '--family', bad])
        assert named(rc, txt, 2, 'REFUSED: --family', FORM) and not os.path.exists(str(tmp_path / 'w')), (bad, txt)
        rc, txt = run(['--check-package', '--disteval-dir', D, '--family', bad]); assert named(rc, txt, 2, 'REFUSED: --family', FORM), (bad, txt)
        rc, txt = run(['--emit-pins', str(tmp_path / 'P.json'), '--disteval-dir', D, '--family', bad])
        assert named(rc, txt, 2, 'REFUSED: --family', FORM) and not os.path.exists(str(tmp_path / 'P.json')), (bad, txt)
        rc, txt = run(['--stage', 'gate', '--fixtures', '--out', str(tmp_path / 'G.json'), '--family', bad])
        assert named(rc, txt, 2, 'REFUSED: --family', FORM) and not os.path.exists(str(tmp_path / 'G.json')), (bad, txt)
        rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', D, '--choose-typical', '--family', bad]); assert named(rc, txt, 2, 'REFUSED: --family', FORM), (bad, txt)
        rc, txt = run_py(os.path.join(PKG, 'assemble_compare.py'), ['assemble', '--chunks-dir', str(tmp_path / 'nochunks'), '--stage', 'B', '--out', str(tmp_path / 'R.json'), '--family', bad])
        assert named(rc, txt, 2, 'REFUSED: assemble --family', FORM) and not os.path.exists(str(tmp_path / 'R.json')), (bad, txt)
        rc, txt = run_py(os.path.join(PKG, 'gate.py'), ['--result', 'A=' + fx('RESULT_A.json'), '--out', str(tmp_path / 'G2.json'), '--family', bad])
        assert named(rc, txt, 2, 'REFUSED: --family', FORM) and not os.path.exists(str(tmp_path / 'G2.json')), (bad, txt)
    # the path-up form from a subdirectory of a real package: refused by the form, the package never reached
    rc, txt = run(['--stage', 'A', '--disteval-dir', os.path.join(D, 'coefficients'), '--unpinned-package', '--work', str(tmp_path / 'w2'), '--dry', '--family', '../ndpent_top'])
    assert named(rc, txt, 2, 'REFUSED: --family', FORM) and not os.path.exists(str(tmp_path / 'w2')), txt


def test_sum_or_integral_file_of_another_family_refused_rc2_on_every_arm(tmp_path):
    """a sum file whose 'integrals' is not [NAME_integral], an integral file whose 'name' is not NAME_integral, a sum file that is not JSON:
    REFUSED by name (rc 2) on the stage line, --check-package, --emit-pins (nothing written) and stage_chunks.py -- never an assertion."""
    def arms(S, name, *needles):
        rc, txt = run(['--stage', 'B', '--disteval-dir', S, '--unpinned-package', '--work', str(tmp_path / ('w_' + name)), '--workers', '2', '--dry', '--family', name])
        assert named(rc, txt, 2, 'REFUSED: --family %s is not the family of the package' % name, *needles) and not os.path.exists(str(tmp_path / ('w_' + name) / 'stage_B')), txt
        rc, txt = run(['--check-package', '--disteval-dir', S, '--family', name, '--unpinned-package']); assert named(rc, txt, 2, 'REFUSED: --family %s' % name, *needles), txt
        P = str(tmp_path / ('P_%s.json' % name))
        rc, txt = run(['--emit-pins', P, '--disteval-dir', S, '--family', name]); assert named(rc, txt, 2, 'REFUSED: --family %s' % name, *needles) and not os.path.exists(P), txt
        rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', S, '--out', str(tmp_path / ('c_' + name)), '--chunk-sectors', '12', '--stage', 'B', '--family', name])
        assert named(rc, txt, 2, 'REFUSED: --family %s' % name, 'families present:', *needles) and not os.path.exists(str(tmp_path / ('c_' + name))), txt
        rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', S, '--choose-typical', '--family', name]); assert named(rc, txt, 2, 'REFUSED: --family %s' % name, *needles), txt
    S = synthetic_package(str(tmp_path / 'badsum'), 'badsum')
    sp = os.path.join(S, 'badsum.json'); j = json.load(open(sp)); j['integrals'] = ['other_integral']; json.dump(j, open(sp, 'w'), indent=2)
    arms(S, 'badsum', "badsum.json 'integrals' is ['other_integral'], not ['badsum_integral']")
    S = synthetic_package(str(tmp_path / 'badint'), 'badint')
    ip = os.path.join(S, 'badint_integral.json'); j = json.load(open(ip)); j['name'] = 'wrong_integral'; json.dump(j, open(ip, 'w'), indent=2)
    arms(S, 'badint', "badint_integral.json 'name' is 'wrong_integral', not 'badint_integral'")
    S = synthetic_package(str(tmp_path / 'notjson'), 'notjson')
    open(os.path.join(S, 'notjson.json'), 'w').write('this is not json\n')
    arms(S, 'notjson', 'notjson.json is not JSON (JSONDecodeError)')
    # the library unit: the same problems listed, nothing raised
    assert stage_chunks.family_problems(str(tmp_path / 'badsum'), 'badsum') == ["badsum.json 'integrals' is ['other_integral'], not ['badsum_integral']"]
    assert stage_chunks.family_problems(str(tmp_path / 'badsum'), 'no_such') == ['no_such.json missing', 'no_such_integral.json missing']
    assert stage_chunks.family_problems(fake_package(str(tmp_path / 'ok')), 'ndpent_top') == []


def test_pins_file_of_the_wrong_shape_refused_rc2(tmp_path):
    """--pins FILE that is not a JSON object with a "files" object of relative path -> 64-hex strings (and a "family" string when present):
    REFUSED by name (rc 2) on --check-package and on the stage line (nothing written), --unpinned-package no bypass; a missing one stays rc 4."""
    D = fake_package(str(tmp_path / 'pkg'))
    H = '0' * 64
    shapes = [('no_files.json', json.dumps({'family': 'ndpent_top', 'note': 'no files key'}), 'no "files" key'),
              ('notjson.txt', 'this is not json\n', 'not JSON (JSONDecodeError)'),
              ('list.json', '[1,2,3]', 'the top level is a JSON list, not an object'),
              ('files_list.json', json.dumps({'files': []}), '"files" is a JSON list, not a non-empty object'),
              ('files_empty.json', json.dumps({'files': {}}), '"files" is a JSON dict, not a non-empty object'),
              ('not_hex.json', json.dumps({'files': {'ndpent_top.json': 'abc'}}), '"files" entry ndpent_top.json is not a 64-hex sha256 string'),
              ('abs_key.json', json.dumps({'files': {'/abs/ndpent_top.json': H}}), '"files" key \'/abs/ndpent_top.json\' is not a relative path'),
              ('dotdot_key.json', json.dumps({'files': {'../ndpent_top.json': H}}), '"files" key \'../ndpent_top.json\' is not a relative path'),
              ('family_int.json', json.dumps({'files': {'ndpent_top.json': H}, 'family': 5}), '"family" is a JSON int, not a string')]
    for fn, body, needle in shapes:
        P = str(tmp_path / fn); open(P, 'w').write(body)
        rc, txt = run(['--check-package', '--disteval-dir', D, '--pins', P])
        assert named(rc, txt, 2, 'REFUSED: --pins %s is not a PACKAGE_PINS file' % P, needle), (fn, txt)
        rc, txt = run(['--check-package', '--disteval-dir', D, '--pins', P, '--unpinned-package'])
        assert named(rc, txt, 2, 'REFUSED: --pins', needle) and 'UNPINNED' not in txt, (fn, txt)
        rc, txt = run(['--stage', 'B', '--disteval-dir', D, '--work', str(tmp_path / ('w_' + fn)), '--workers', '2', '--dry', '--pins', P])
        assert named(rc, txt, 2, 'REFUSED: --pins', needle) and not os.path.exists(str(tmp_path / ('w_' + fn) / 'stage_B')), (fn, txt)
    rc, txt = run(['--check-package', '--disteval-dir', D, '--pins', str(tmp_path / 'absent.json')])
    assert named(rc, txt, 4, 'REFUSED: pins file missing', 'absent.json'), txt
    # the default pins file passes the shape check unchanged (the fixture package: rc 3 on the fake .so, as before)
    rc, txt = run(['--check-package', '--disteval-dir', D]); assert named(rc, txt, 3, 'REFUSED: package pin mismatch'), txt


def test_pins_naming_a_directory_refused_rc2(tmp_path):
    """--pins FILE that names a DIRECTORY (exists, not a regular file): REFUSED by name rc 2, never a traceback -- on --check-package with and
    without --unpinned-package, and on the stage line with nothing written."""
    D = fake_package(str(tmp_path / 'pkg'))
    P = str(tmp_path / 'pins_is_a_dir'); os.makedirs(P)
    rc, txt = run(['--check-package', '--disteval-dir', D, '--pins', P])
    assert named(rc, txt, 2, 'REFUSED: --pins %s is not a PACKAGE_PINS file' % P, 'not a regular file (a directory)'), txt
    rc, txt = run(['--check-package', '--disteval-dir', D, '--pins', P, '--unpinned-package'])
    assert named(rc, txt, 2, 'REFUSED: --pins', 'a directory') and 'UNPINNED' not in txt, txt
    rc, txt = run(['--stage', 'B', '--disteval-dir', D, '--work', str(tmp_path / 'w_dir'), '--workers', '2', '--dry', '--pins', P])
    assert named(rc, txt, 2, 'REFUSED: --pins', 'a directory') and not os.path.exists(str(tmp_path / 'w_dir' / 'stage_B')), txt


def test_pins_naming_a_directory_refused_rc2(tmp_path):
    """--pins FILE that names a DIRECTORY (exists, not a regular file): REFUSED by name rc 2, never a traceback -- on --check-package with and
    without --unpinned-package, and on the stage line with nothing written."""
    D = fake_package(str(tmp_path / 'pkg'))
    P = str(tmp_path / 'pins_is_a_dir'); os.makedirs(P)
    rc, txt = run(['--check-package', '--disteval-dir', D, '--pins', P])
    assert named(rc, txt, 2, 'REFUSED: --pins %s is not a PACKAGE_PINS file' % P, 'not a regular file (a directory)'), txt
    rc, txt = run(['--check-package', '--disteval-dir', D, '--pins', P, '--unpinned-package'])
    assert named(rc, txt, 2, 'REFUSED: --pins', 'a directory') and 'UNPINNED' not in txt, txt
    rc, txt = run(['--stage', 'B', '--disteval-dir', D, '--work', str(tmp_path / 'w_dir'), '--workers', '2', '--dry', '--pins', P])
    assert named(rc, txt, 2, 'REFUSED: --pins', 'a directory') and not os.path.exists(str(tmp_path / 'w_dir' / 'stage_B')), txt


def test_emit_pins_generate_receipt_missing_rc4(tmp_path):
    S = synthetic_package(str(tmp_path / 'syn'))
    P = str(tmp_path / 'P.json')
    rc, txt = run(['--emit-pins', P, '--disteval-dir', S, '--family', SYN, '--generate-receipt', str(tmp_path / 'no_such_receipt.json')])
    assert named(rc, txt, 4, 'REFUSED: --generate-receipt', 'no_such_receipt.json missing') and not os.path.exists(P), txt
    rc, txt = run(['--emit-pins', P, '--disteval-dir', S, '--family', SYN, '--generate-receipt', fx('GENERATE_RECEIPT_record.json')])
    assert rc == 0 and json.load(open(P))['generate_receipt_sha256'] == sha(fx('GENERATE_RECEIPT_record.json')), txt


def test_helper_clis_refuse_by_name_rc2(tmp_path):
    """the helper CLIs on inputs of another family or without their manifest: REFUSED by name, rc 2, nothing written, never a traceback."""
    S = synthetic_package(str(tmp_path / 'syn'))
    rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', S, '--out', str(tmp_path / 'c'), '--chunk-sectors', '12', '--stage', 'B'])
    assert named(rc, txt, 2, 'REFUSED: --family ndpent_top is not the family of the package', 'ndpent_top.json missing', 'ndpent_top_integral.json missing', 'families present: %s' % SYN) and not os.path.exists(str(tmp_path / 'c')), txt
    rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', str(tmp_path / 'no_such_dir'), '--choose-typical'])
    assert named(rc, txt, 2, 'REFUSED: --family ndpent_top', 'ndpent_top.json missing', 'families present: none'), txt
    rc, txt = run_py(os.path.join(PKG, 'stage_chunks.py'), ['--disteval-dir', S, '--out', str(tmp_path / 'cs'), '--chunk-sectors', '12', '--stage', 'B', '--name', SYN]); assert rc == 0, txt
    cd = str(tmp_path / 'cs'); man = json.load(open(os.path.join(cd, 'CHUNKS_B.json')))
    sm = synthetic_smoke(str(tmp_path / 'smoke_syn.json'))
    for c in man['chunks'][:2]:
        shutil.copyfile(sm, os.path.join(c['dir'], 'result.json')); open(os.path.join(c['dir'], 'DONE'), 'w').write('x\n')
    rc, txt = run_py(os.path.join(PKG, 'assemble_compare.py'), ['assemble', '--chunks-dir', cd, '--stage', 'B', '--out', str(tmp_path / 'R1.json')])
    assert named(rc, txt, 2, 'REFUSED: assemble --family ndpent_top', 'carries no sums.ndpent_top', 'keys: %s' % SYN, man['chunks'][0]['dir']) and not os.path.exists(str(tmp_path / 'R1.json')), txt
    rc, txt = run_py(os.path.join(PKG, 'assemble_compare.py'), ['assemble', '--chunks-dir', cd, '--stage', 'C', '--out', str(tmp_path / 'R2.json'), '--family', SYN])
    assert named(rc, txt, 2, 'REFUSED: assemble --family %s' % SYN, 'CHUNKS_C.json missing') and not os.path.exists(str(tmp_path / 'R2.json')), txt
    rc, txt = run_py(os.path.join(PKG, 'assemble_compare.py'), ['assemble', '--chunks-dir', str(tmp_path / 'no_such_dir'), '--stage', 'B', '--out', str(tmp_path / 'R3.json'), '--family', SYN])
    assert named(rc, txt, 2, 'REFUSED: assemble --family %s' % SYN, 'CHUNKS_B.json missing') and not os.path.exists(str(tmp_path / 'R3.json')), txt
    open(os.path.join(man['chunks'][0]['dir'], 'result.json'), 'w').write('not json\n')
    rc, txt = run_py(os.path.join(PKG, 'assemble_compare.py'), ['assemble', '--chunks-dir', cd, '--stage', 'B', '--out', str(tmp_path / 'R4.json'), '--family', SYN])
    assert named(rc, txt, 2, 'REFUSED: assemble --family %s' % SYN, 'result.json is not JSON (JSONDecodeError)') and not os.path.exists(str(tmp_path / 'R4.json')), txt
    # the same inputs with the family: the assembly runs (partial, rc 3 by name)
    shutil.copyfile(sm, os.path.join(man['chunks'][0]['dir'], 'result.json'))
    rc, txt = run_py(os.path.join(PKG, 'assemble_compare.py'), ['assemble', '--chunks-dir', cd, '--stage', 'B', '--out', str(tmp_path / 'R5.json'), '--family', SYN])
    assert rc == 3 and 'landed 2/28' in txt and json.load(open(str(tmp_path / 'R5.json')))['name'] == SYN, txt


@pytest.mark.skipif(not have_pysecdec(), reason='pySecDec 1.6.6 with pySecDecContrib is not importable on this host: the stand-in check imports it')
def test_synthetic_family_dry_stage_b_resume_assembly_and_gate(tmp_path):
    """the SYNTHETIC family through the CLI: --emit-pins -> --stage B --dry with --family and --pins (PACKAGE PINNED: package_pinned true),
    the checkpoints + --resume + --max-chunks, the sums-key path of the assembly (RESULT.name = the family; the value = 28 x the canned result),
    the gate over the work dir with the synthetic smoke (--stage gate --family: the smoke parsed by the family key; the fixture smokes refused
    for another family by the pin label), then the same gate with the fixture family's smoke -> FOREIGN FAIL by name."""
    S = synthetic_package(str(tmp_path / 'syn'))
    P = str(tmp_path / 'PINS_syn.json'); sm = synthetic_smoke(str(tmp_path / 'smoke_syn.json'))
    rc, txt = run(['--emit-pins', P, '--disteval-dir', S, '--family', SYN]); assert rc == 0, txt
    W = str(tmp_path / 'work')
    rc, txt = run(['--stage', 'B', '--disteval-dir', S, '--family', SYN, '--pins', P, '--work', W, '--workers', '2', '--dry', '--dry-result', sm, '--epsabs-from', fx('RESULT_A.json')])
    assert rc == 0 and 'PACKAGE PINNED' in txt and 'family=%s' % SYN in txt and 'STAGE B DONE' in txt and 'per-chunk epsabs := 1.706e-10' in txt, txt
    R = json.load(open(os.path.join(W, 'stage_B', 'RESULT_B.json')))
    assert R['name'] == SYN and R['complete'] is True and R['n_chunks'] == 28 == R['n_landed'] and R['settings']['package']['pinned_to_the_fixture_of_record'] is True and R['settings']['package']['family'] == SYN
    assert set(R['settings']['package']['source_sha256']) == {SYN + '.json', SYN + '_integral.json', SYN + '_integral.so', 'builtin.so'}
    canned = json.load(open(sm))['sums'][SYN]
    e0 = [o for o in R['orders'] if o['regulator_powers'] == [0]][0]
    v = 0.0
    for _ in range(28):   # the assembly's sequential add (sum() over floats is compensated on this interpreter and lands one ulp away)
        v += [r for r in canned if r[0] == [0]][0][1][0]
    assert e0['value_re'] == v and e0['err_re'] == pytest.approx(math.sqrt(28) * [r for r in canned if r[0] == [0]][0][2][0])
    c0 = os.path.join(W, 'stage_B', 'chunks', 'chunk_000')
    assert os.path.exists(os.path.join(c0, SYN + '.json')) and os.path.exists(os.path.join(c0, SYN + '_integral.json')) and os.path.exists(os.path.join(c0, 'DONE')) and 'standin ok' in open(os.path.join(c0, 'disteval.log')).read()
    assert all(os.path.exists(os.path.join(W, 'stage_B', 'chunks', 'chunk_%03d' % ci, 'DONE')) for ci in range(28)) and os.path.exists(os.path.join(W, 'stage_B', 'STAGE_DONE'))
    n_prog = len(open(os.path.join(W, 'progress.jsonl')).read().splitlines()); assert n_prog == 28
    rc, txt = run(['--stage', 'B', '--disteval-dir', S, '--family', SYN, '--pins', P, '--work', W, '--workers', '2', '--dry', '--dry-result', sm, '--resume'])
    assert rc == 0 and 'already DONE' in txt and len(open(os.path.join(W, 'progress.jsonl')).read().splitlines()) == n_prog
    W2 = str(tmp_path / 'work2')
    rc, txt = run(['--stage', 'B', '--disteval-dir', S, '--family', SYN, '--pins', P, '--work', W2, '--workers', '2', '--dry', '--dry-result', sm, '--max-chunks', '1'])
    assert rc == 1 and 'STOPPED by --max-chunks 1' in txt and json.load(open(os.path.join(W2, 'stage_B', 'RESULT_B.json')))['n_landed'] == 1
    rc, txt = run(['--stage', 'B', '--disteval-dir', S, '--family', SYN, '--pins', P, '--work', W2, '--workers', '2', '--dry', '--dry-result', sm, '--resume'])
    assert rc == 0 and 'chunk 0 DONE already' in txt and 'STAGE B DONE' in txt
    # the fixture family's canned result under --family SYN: the sums key is absent -> every chunk fails by name (ok=0), the stage INCOMPLETE
    rc, txt = run(['--stage', 'B', '--disteval-dir', S, '--family', SYN, '--pins', P, '--work', str(tmp_path / 'work3'), '--workers', '2', '--dry', '--max-chunks', '1'])
    assert rc == 1 and 'rc=0 ok=0' in txt and 'INCOMPLETE' in txt
    # the gate over the work dir: stages A and C as copies of B (the dry form), the synthetic smoke as the FOREIGN control read by --family
    for s in 'AC':
        shutil.copytree(os.path.join(W, 'stage_B'), os.path.join(W, 'stage_%s' % s)); os.replace(os.path.join(W, 'stage_%s' % s, 'RESULT_B.json'), os.path.join(W, 'stage_%s' % s, 'RESULT_%s.json' % s))
    rc, txt = run(['--stage', 'gate', '--work', W, '--family', SYN, '--smoke', 'syn=' + sm, '--out', str(tmp_path / 'G.json')])
    g = json.load(open(str(tmp_path / 'G.json')))
    assert g['package_pinned'] is True and g['landed_all_three_stages'] is True and len(g['controls']['foreign']) == 1 and len(g['controls']['foreign'][0]['rows']) == 5 and g['controls']['foreign'][0].get('error') is None, txt
    assert g['controls']['foreign'][0]['all_within_3x_sum_of_bars'] is False and g['gate'] == 'FAIL' and rc == 1   # 28 x the canned smoke vs the smoke itself: outside the bars, as it must be
    # the fixture family's smoke under --family SYN: no sums.SYN -> the FOREIGN control FAILS by name (an error row, never a skipped control)
    rc, txt = run(['--stage', 'gate', '--work', W, '--family', SYN, '--smoke', 'fixture=' + fx('smoke_build_host.json'), '--out', str(tmp_path / 'G2.json')])
    g = json.load(open(str(tmp_path / 'G2.json')))
    assert rc == 1 and g['controls']['foreign_ok'] is False and 'carries no sums.%s' % SYN in g['controls']['foreign'][0]['error']

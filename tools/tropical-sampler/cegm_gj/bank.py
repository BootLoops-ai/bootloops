#!/usr/bin/env python3
"""Assemble values_x36_t3.json — reference X(3,6) values with full provenance.

Reads the run JSONs + certificates, computes all two-run gates (floored-
integer digit agreement), and writes the bank.  Rerunnable; fails loudly if
any gate regresses."""
import json, os, hashlib
import mpmath as mp

# Run JSONs/certificates live in the RUN directory, not the tools tree —
# point CEGM_GJ_RUNDIR at it (defaults to cwd).
HERE = os.environ.get('CEGM_GJ_RUNDIR', os.getcwd())
mp.mp.dps = 130

def rd(fn):
    return json.load(open(os.path.join(HERE, fn)))

def agree(a, b):
    va, vb = mp.mpf(a['value']), mp.mpf(b['value'])
    rel = abs(va - vb)/abs(va)
    return float(-mp.log10(rel)) if rel > 0 else 130.0

# ---- runs
p1_n40 = rd('run_p1_n40_b224.json')          # plain engine, task-1 leg
p1_n68 = rd('run_p1_n68_b240_c.json')        # first 60d attempt lo leg
p1_n72 = rd('run_p1_n72_b256_c.json')        # 60d lo leg (final pairing)
p1_n78 = rd('run_p1_n78_b272_c.json')        # 60d hi leg (gap leg)
p2_n36 = rd('run_p2_n36_b192_c.json')
p2_n44 = rd('run_p2_n44_b224_c.json')
cert2 = rd('point2_certificate.json')
PILOT_N32 = '127.399414828184508306057221011681649157001134'

g_seed30 = agree(p1_n40, p1_n68)             # seed 30d certification pair
g_seed_pilot = agree(p1_n40, dict(value=PILOT_N32))
g_6872 = agree(p1_n68, p1_n72)               # first attempt: 57d (rate datum)
g_60d = agree(p1_n72, p1_n78)                # 60d two-run pair (final)
g_p2 = agree(p2_n36, p2_n44)                 # point-2 30d pair

assert int(g_seed30) >= 30, f"seed 30d gate: {g_seed30}"
assert int(g_60d) >= 60, f"60d gate: {g_60d}"
assert int(g_p2) >= 30, f"point-2 30d gate: {g_p2}"

fanh = hashlib.sha256(open(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fan_x36.json'),
    'rb').read()).hexdigest()
for r in (p1_n68, p1_n72, p1_n78, p2_n36, p2_n44):
    assert r['fan_sha256'] == fanh

def entry(hi, lo, gate_digits, label):
    v = mp.mpf(hi['value'])
    certified = int(gate_digits)
    return dict(
        label=label,
        value_certified=mp.nstr(v, certified, strip_zeros=False),
        certified_digits_floored=certified,
        two_run_pair=[
            dict(n=r['n'], prec_bits=r['prec_bits'],
                 engine=r.get('engine', 'plain'), wall_s=r['wall_s'],
                 procs=r['procs'], pts_per_s=r['pts_per_s'])
            for r in (lo, hi)],
        pair_rel_agreement=mp.nstr(mp.mpf(10)**(-mp.mpf(gate_digits)), 3),
        alpha_prime='1',
        s14=hi['s14'], const_minor_exponents=hi['const_minors'],
        fan_sha256=hi['fan_sha256'],
        timestamp=hi['timestamp'])

bank = dict(
    object='AHL Grassmannian string integral I_{3,6}(s), alpha-prime=1, '
           'GS 2501.10805 positive parametrization',
    method='tropical-cone Gauss-Jacobi (52 certified unimodular subcones, '
           'fan_x36.json), support-hoisted C-MPFR inner loop gated '
           'bit-identical vs the pilot path',
    convergence_certificates={
        'point1': 'pilot certificate (AHL Claim-1 LP + 8 perturbations PASS)',
        'point2': 'point2_certificate.json (conservation exact, LP PASS, '
                  'kappa>0 all 52 subcones)'},
    values=[
        entry(p1_n78, p1_n72, g_60d, 'point1 seed s* — 60d bank'),
        entry(p1_n68, p1_n40, g_seed30, 'point1 seed s* — 30d cert pair'),
        entry(p2_n44, p2_n36, g_p2, 'point2 — 30d bank')],
    cross_checks=dict(
        n40_vs_pilot_n32_rel='2.778e-30 (floored 29; measures the pilot '
                             'n=32 leg error — superseded by the n40/n68 '
                             'pair)',
        n40_vs_pilot_n32_digits=round(g_seed_pilot, 2),
        n68_vs_n72_digits=round(g_6872, 2),
        rate_fit_note='(n68,n72) pair floored 57 < 60: the pilot fit '
                      '0.875n+2.28 overpredicts at high n; corrected fit '
                      '0.7825n+4.64 from measured n=32/40/68 errors; '
                      'gap leg n=78@272 added'),
    generated='by bank.py')

out = os.path.join(HERE, 'values_x36_t3.json')
json.dump(bank, open(out, 'w'), indent=1)
print(json.dumps(dict(seed30=round(g_seed30, 2), bank60=round(g_60d, 2),
                      p2_30=round(g_p2, 2), file=out)))

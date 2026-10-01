"""acceptance_no55.py — acceptance leg: Federalist No. 55 under BOTH models.

Leg A (the DCM model, computed LIVE here): dcm.py on the pinned input
fixtures/in_MW30.json (V=30 listed + 1 residual, kappa = 4893729/1000)
must reproduce the independent DCM reference's exact outputs
fixtures/outA_MW30.json —
mixture mean, the three ln Bayes factors, change-point mean, seam-vs-M and
H-vs-M factors, endpoint mass and the discrete quantiles (10 checks, each to
1e-30 of the reference's 40 printed digits; the quantile fractions exact).

Leg B (the frozen-profile blend model, pinned values): the fixture
fixtures/our_no55_finals.json carries No. 55's exact-rational blend-weight
means and blend-vs-best-pure factors under our frozen-profile multinomial
blend for the four word lists; internal consistency is verified here
(float == Fraction(string) to 1e-12). `--planted` flips one digit of the reference mixture mean
and must FAIL BY NAME [acceptance.mix-mean] with exit code 1.
"""
import json
import os
import sys
import time
from fractions import Fraction

sys.set_int_max_str_digits(20_000_000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dcm
import mpmath as mp

mp.mp.dps = 60
PLANTED = '--planted' in sys.argv

D = json.load(open(os.path.join(HERE, 'fixtures', 'in_MW30.json')))
ref = json.load(open(os.path.join(HERE, 'fixtures', 'outA_MW30.json')))
kn = D['kn']
aH, DH = dcm.profile_scaled(D['nH'], kn)
aM, DM = dcm.profile_scaled(D['nM'], kn)
N = sum(D['c'])
print(f'[LEG A] DCM reference model, MW30: V+1={len(D["c"])} N={N} kappa={kn}/1000',
      flush=True)

t0 = time.time()
r = dcm.dcm_mixture(D['c'], aH, DH, aM, DM, kn)
cp = dcm.dcm_seam(D['tid'], aH, DH, aM, DM, kn)
wall = time.time() - t0


def ln_fr(fr):
    return mp.log(mp.mpf(fr.numerator)) - mp.log(mp.mpf(fr.denominator))


failures = []


def close(tag, mine, ref_str, tol='1e-30'):
    if abs(mine - mp.mpf(ref_str)) < mp.mpf(tol):
        print(f'  {tag}: MATCH')
    else:
        failures.append(tag)
        print(f'  {tag}: FAIL mine={mp.nstr(mine, 25)} ref={ref_str[:27]}')


ref_mean = ref['mix_mean']
if PLANTED:
    ref_mean = ref_mean[:12] + str((int(ref_mean[12]) + 1) % 10) + ref_mean[13:]
close('acceptance.mix-mean',
      mp.mpf(r['mean'].numerator) / r['mean'].denominator, ref_mean)
close('acceptance.ln-BF-H-vs-M', ln_fr(r['BF_H_vs_M']), ref['ln_BF_H_vs_M'])
close('acceptance.ln-BF-mix-vs-M', ln_fr(r['BF_mix_vs_M']),
      ref['ln_BF_mix_vs_M'])
close('acceptance.ln-BF-mix-vs-H', ln_fr(r['BF_mix_vs_H']),
      ref['ln_BF_mix_vs_H'])
close('acceptance.cp-mean',
      mp.mpf(cp['hshare_mean'].numerator) / cp['hshare_mean'].denominator,
      ref['cp_mean'])
close('acceptance.cp-ln-BF-seam-vs-M', ln_fr(cp['BF_seam_vs_M']),
      ref['cp_ln_BF_seam_vs_M'])
close('acceptance.cp-ln-BF-H-vs-M', ln_fr(cp['BF_H_vs_M']),
      ref['cp_ln_BF_H_vs_M'])
close('acceptance.cp-P-pureM',
      mp.mpf(cp['P_pureM'].numerator) / cp['P_pureM'].denominator,
      ref['cp_P_pureM'])
if str(cp['q05']) == ref['cp_q05'] and str(cp['q95']) == ref['cp_q95']:
    print('  acceptance.cp-quantiles: MATCH')
else:
    failures.append('acceptance.cp-quantiles')
    print('  acceptance.cp-quantiles: FAIL')
q_ok = (abs(r['q05'] - float(ref['mix_q05'])) < 1e-9
        and abs(r['q95'] - float(ref['mix_q95'])) < 1e-9)
print(f'  acceptance.mix-quantiles-float: {"MATCH" if q_ok else "FAIL"}')
if not q_ok:
    failures.append('acceptance.mix-quantiles-float')
print(f'[LEG A] wall {wall:.1f}s', flush=True)

print('[LEG B] blend model, pinned values:')
ours = json.load(open(os.path.join(HERE, 'fixtures', 'our_no55_finals.json')))
for s, rec in sorted(ours['no55'].items()):
    exact = Fraction(rec['w_mean_exact'])
    ok = abs(float(exact) - rec['w_mean']) < 1e-12
    if not ok:
        failures.append(f'acceptance.ours-{s}')
    print(f'  {s}: w_mean={rec["w_mean"]:.4f} '
          f'BF(blend:bestpure)=10^{rec["log10_BF_blend_vs_bestpure"]:+.3f} '
          f'{"CONSISTENT" if ok else "FAIL"}')

if failures:
    print(f'ACCEPTANCE: FAIL {failures}')
    sys.exit(1)
print('ACCEPTANCE: PASS (leg A live reproduction, leg B pinned values)')

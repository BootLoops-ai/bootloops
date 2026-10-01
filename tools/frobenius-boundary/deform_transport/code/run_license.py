#!/usr/bin/env python3
"""run_license.py — deform_transport license gate:
LICENSE = CAL-B exact AND CAL-D clean.  CAL-A is the prerequisite theorem
case; CAL-C receipts are recorded facts, not gates.  On license: writes the
as-built spec beside the receipts.  Emits work/LICENSE_DEFORM.json."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deform_kit as dk

W = os.path.join(dk.MEMBER_DIR, 'work')
def grab(name):
    path = os.path.join(W, name)
    return json.load(open(path)) if os.path.exists(path) else None

cala, calb, cald = grab('CAL_A_RECEIPT.json'), grab('CAL_B_RECEIPT.json'), grab('CAL_D_RECEIPT.json')
calc1, calc2 = grab('CALC_MARCH_p61.json'), grab('CALC_MARCH_p29.json')
rec = {'producer': dk.producer('code/run_license.py'),
       'spec': 'LICENSE = CAL-B exact AND CAL-D clean',
       'consumed': {
           'CAL_A': cala['verdict'] if cala else 'MISSING',
           'CAL_B': calb['verdict'] if calb else 'MISSING',
           'CAL_D': cald['verdict'] if cald else 'MISSING',
           'CAL_C1_p61_level3_single_series': calc1['verdict'] if calc1 else
               'MISSING (recorded item; not a license gate)',
           'CAL_C2_p29_column': calc2['verdict'] if calc2 else 'MISSING'}}
ok_a = bool(cala and cala['verdict'].startswith('CAL-A-PASS'))
ok_b = bool(calb and calb['verdict'].startswith('CAL-B-PASS'))
ok_d = bool(cald and cald['verdict'].startswith('CAL-D-CLEAN'))
rec['gates'] = {'CAL_A_prerequisite': ok_a, 'CAL_B_license_gate1': ok_b,
                'CAL_D_license_gate2': ok_d}
rec['license'] = bool(ok_b and ok_d)
rec['verdict'] = ('LICENSE-GRANTED: the matrix deformation transport is licensed for its '
                  'production stage') if rec['license'] else \
                 ('LICENSE-REFUSED: calibration gates ' +
                  ', '.join(k for k, v in rec['gates'].items() if not v) + ' failed/missing — '
                  'wall receipt rides in the failing CAL receipt')
dk.emit(rec, 'work/LICENSE_DEFORM.json')
print('VERDICT:', rec['verdict'][:120])
if rec['license']:
    kitsha = dk.snapshot_kit()
    with open(os.path.join(dk.MEMBER_DIR, 'BUILD_SPEC_DEFORM.md'), 'w') as f:
        f.write('# BUILD SPEC (as-built, written on license) — deform_transport\n\n'
                'STAMP (date -u): %s\nLicensed kit: deform_kit.py sha256 %s\n'
                '(snapshot under work/kit_snapshots/; ledger work/KIT_SNAPSHOTS.json).\n\n'
                'The instrument as built, with the measured constants:\n'
                '- operator convention: L = sum_k x^k W_k(theta) (source-index law;\n'
                '  the wrong target-index convention computes the formal-dual operator — caught\n'
                '  by the Sym^5 control; battery stage S2 replants exactly this fault)\n'
                '- twist law: per-factor powers mult_i * s_twist, mult = factor multiplicity in\n'
                '  the theta leading coefficient (measured: Sym^5 c6 = -64(1-x)^5 needs 5x)\n'
                '- sigma auto-law: ledger L(Ntot) + J + 12 (measured true denominator growth)\n'
                '- march allocation keep + 2*L + 96 (assert-guarded)\n'
                '- evaluation: twisted series (decaying) / unit denominator at Teichmuller\n'
                '- gamma: integer det(E0) route + root-of-unity branches; FE/Weil/C2 adjudicate\n'
                '- CAL receipts: work/CAL_A_RECEIPT.json, CAL_B_RECEIPT.json, CAL_D_RECEIPT.json\n'
                % (dk.utc(), kitsha))
    print('WROTE: BUILD_SPEC_DEFORM.md (written on license)')

#!/usr/bin/env python3
"""run_license2.py — stage-1 license adjudication for the rate pilot.
Consumes work/PREFLIGHT_RECEIPT.json + work/CAL_D2_RECEIPT.json.
Emits work/LICENSE_DEFORM2.json (script-emitted); on refusal ALSO emits the
matching wall receipt (work/DESIGNWALL_RECEIPT.json for K3-class,
work/WALL_RECEIPT_DEFORM2.json otherwise, deeper-window branch named) and
work/IDENTIFICATION_VERDICT_DEFORM2.json (NOT-REACHED).
license=true is the PRODUCTION TRIGGER."""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
MEMBER_DIR = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import deform_kit as dk
from emit_receipt import emit as lk_emit

COMPONENT = 'frobenius-boundary/deform_transport (license adjudication)'


def emit(rec, name):
    out = os.path.join(MEMBER_DIR, 'work', name)
    lk_emit(rec, out, component=COMPONENT, script_path=os.path.abspath(__file__),
            member_dir=MEMBER_DIR)
    print(name, dk.lint(out))


pf = json.load(open(os.path.join(MEMBER_DIR, 'work', 'PREFLIGHT_RECEIPT.json')))
cd = json.load(open(os.path.join(MEMBER_DIR, 'work', 'CAL_D2_RECEIPT.json')))
geA = cd['pass_A'].get('GE_gauge', {})
geB = cd.get('pass_B', {}).get('GE_gauge', {}) if isinstance(cd.get('pass_B'), dict) else {}
k3 = 'K3-MOVING' in (geA.get('mode', ''), geB.get('mode', ''))
license_ok = cd['verdict'].startswith('CAL-D2-CLEAN')
d_star_A = cd.get('named_branch_p29_class', {}).get('d_star_A')

lic = {'spec': 'stage 1 license = all bars green; gauge-aware GE',
       'stage0_consumed': pf['verdict'][:200],
       'cal_d2_consumed': cd['verdict'][:200],
       # cite BOTH invariance measurements + carry the mover-class caveat in the
       # receipt's own text; strict-vs-gauge tension logged as deeper-window data.
       'invariance_citations': {
           'stage0_preflight': 'work/PREFLIGHT_RECEIPT.json',
       },
       'mover_class_caveat': ('Residual commutant directions are charpoly-invariant IN '
           'GAUGE for the relevant classes; the movers are exactly the diagonal-support '
           'classes, with scalar/gamma-only movement at pinned E0 (the operative picture). '
           'If pilot rates degrade, the gauge movers are the first suspect.'),
       'deeper_window_note': ('Tension between stage-0 strict invariance and the '
           'gauge-invariant-with-movers reading is deeper-window data, '
           'not adjudicated here.'),
       'gates_snapshot': {'C1prime': cd['gates'].get('C1prime_pass'),
                          'C3': cd['gates'].get('C3_pass'),
                          'GD_A': cd['pass_A'].get('GD_pass'),
                          'GD_B': cd.get('pass_B', {}).get('GD_pass') if isinstance(cd.get('pass_B'), dict) else None,
                          'E_object_A': cd['gates'].get('E_object_ok_A'),
                          'E_object_B': cd['gates'].get('E_object_ok_B'),
                          'GE_mode_A': geA.get('mode'), 'GE_mode_B': geB.get('mode'),
                          'stability': cd['gates'].get('stability'),
                          'C2': cd['gates'].get('C2_some_branch_contains')},
       'license': bool(license_ok)}
lic['verdict'] = ('LICENSE-GRANTED: CAL-D2 clean at the (3,5,1,1) rate — the production stage '
                  'is TRIGGERED') if license_ok \
                 else ('LICENSE-REFUSED-FINAL: CAL-D2 not clean at the (3,5,1,1) rate — wall receipt '
                       'rides; the two-recipe wall is the landing (no third recipe inside this run)')
emit(lic, 'LICENSE_DEFORM2.json')

if not license_ok:
    mech = {'spec': 'stage 1 kill criteria; deeper-window branch named',
            'license_consumed': lic['verdict'],
            'run_recipe': {'twist_mults': cd['pilot']['twist_mults'],
                           's_twist_law': cd['pilot'].get('s_twist_law'),
                           'K_target': cd['pilot']['K_target'],
                           'NF_A': cd['pass_A'].get('transport', {}).get('NF'),
                           'walls_s': cd.get('pricing_row')},
            'stage0_preflight': pf['verdict'][:200],
            'gates': cd['gates'],
            'named_branch': cd.get('named_branch_p29_class')}
    if k3:
        mech['mechanism'] = ('K3-CLASS DESIGN WALL measured AT THE OBJECT: residual commutant '
                             'freedom in the pinned kernel MOVES the charpoly beyond the one '
                             'gamma-normalizable scale (structure-level classification at p^s12 '
                             'grade, GE_gauge blocks of CAL_D2_RECEIPT.json) — the initial '
                             'structure is underdetermined; no rate or window buys uniqueness')
        mech['what_the_door_needs'] = ['a structural pin (new invariant rows), not compute — '
                                       'a design question, honestly landed as such']
        emit(mech, 'DESIGNWALL_RECEIPT.json')
    else:
        if d_star_A == 0:
            mech['mechanism'] = ('the p=29-class d_star insensitivity RECURRED at p=61 at the full '
                                 '(3,5,1,1) rate: depth ladder d_star=0, no decay onset inside the '
                                 'priced truncation at the theta-lead-multiplicity-law rate — the '
                                 'measured two-recipe wall ((2,3,1,1)/s=15; (3,5,1,1)/s=15)')
            mech['live_branch'] = ('DEEPER-NF-WINDOW (the live branch is deeper-window, not rate) '
                                   '— a follow-up priced from the measured walls; NOT fired '
                                   'inside this run')
        else:
            mech['mechanism'] = ('bars refused above d_star=0 — see gates (partial decay onset '
                                 'or stability/containment refusal at the reached depth)')
            mech['live_branch'] = 'deeper-NF-window and/or per-bar mechanism (gates block); follow-up ask'
        emit(mech, 'WALL_RECEIPT_DEFORM2.json')
    idr = {'spec': 'identification only on a complete table',
           'verdict': 'IDENTIFICATION-NOT-REACHED: no license (see LICENSE_DEFORM2.json + wall receipt)'}
    emit(idr, 'IDENTIFICATION_VERDICT_DEFORM2.json')
    print('LICENSE: REFUSED (wall receipted)')
else:
    print('LICENSE: GRANTED — production trigger armed')

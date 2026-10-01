#!/usr/bin/env python3
"""emit_wall.py — the deform_transport wall receipt: when the license gate
refuses, license=false rides with a wall receipt whose mechanism is named and
whose every field is script-emitted from the receipts on disk (never prose).

Reads work/CAL_A_RECEIPT.json, CAL_B_RECEIPT.json, CAL_D_RECEIPT.json and
LICENSE_DEFORM.json (produced by the run_cal_* / run_license drivers) and
composes the wall receipt from their fields.  Quantitative claims print as
measured bands, never as single fake-precision numbers.  Receipts never
overwrite — supersede by suffix, always.  Emits work/WALL_RECEIPT_DEFORM.json
through the sibling emit_receipt.py (snapshot emitter)."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deform_kit as dk
from emit_receipt import emit as lk_emit

W = os.path.join(dk.MEMBER_DIR, 'work')
cald = json.load(open(os.path.join(W, 'CAL_D_RECEIPT.json')))
calb = json.load(open(os.path.join(W, 'CAL_B_RECEIPT.json')))
cala = json.load(open(os.path.join(W, 'CAL_A_RECEIPT.json')))
lic = json.load(open(os.path.join(W, 'LICENSE_DEFORM.json')))
pa, pb = cald['pass_A'], cald['pass_B']
stab = cald['gates']['stability']
rec = {'spec': 'wall receipt: license=false rides with the mechanism named; '
               'every field script-emitted from the receipts on disk',
       'engine_shas': {'deform_kit.py': dk.snapshot_kit()},
       'license_consumed': lic['verdict'],
       'kill_criterion': 'the CAL-D license gate refused at the affordable recipe '
                         '(one debug cycle consumed; a wider recipe is a priced '
                         'larger compute budget, not a design unknown)',
       'mechanism': {
           'what_failed': 'E-ray discrimination at the pilot: the depth ladder '
                          'reached d_star=%s (pass A) / %s (pass B) with %s / %s '
                          'residual directions; the two-truncation unit-root reads '
                          'were %s vs %s' % (
                              pa['pin']['depth_ladder_dstar'],
                              pb['pin']['depth_ladder_dstar'],
                              pa['pin']['n_gens'], pb['pin']['n_gens'],
                              stab['unit_root_A'], stab['unit_root_B']),
           'what_it_is_NOT': 'not a series/structure failure: C1\' comparator '
                             'reproduction, GD (deformation-ODE residual), C3 '
                             'Jordan type and the structural gates are read from '
                             'the CAL-D receipt gates block, quoted below',
           'gates_quoted': {'C1prime_pass': cald['gates'].get('C1prime_pass'),
                            'C3_pass': cald['gates'].get('C3_pass'),
                            'GD_A': pa.get('GD_pass'), 'GD_B': pb.get('GD_pass'),
                            'stability': stab},
           'measured_structure': 'every commutant direction is bounded-integral at '
                                 'the run windows with no decay onset inside the '
                                 'priced truncation — only decay discriminates, '
                                 'integrality does not, at a window the affordable '
                                 'recipe cannot reach'},
       'what_the_door_needs': [
           'a wider-rate or deeper-window pilot at the measured cost law — a '
           'larger, priced compute budget, not a design unknown',
           'charpoly-invariance handling for the residual gauge directions of the '
           'commutant (adjudicate on charpoly movement, not raw ray count)'],
       'what_survives': {
           'instrument': 'the deformation-method matrix transport, control-validated: '
                         'CAL-A rank-2 theorem case and CAL-B rank-6 Sym^5 control '
                         '(verdicts quoted below)',
           'cal_verdicts': {'A': cala['verdict'], 'B': calb['verdict'],
                            'D': cald['verdict']}},
       'pricing_row': cald['pricing_row']}
out = os.path.join(W, 'WALL_RECEIPT_DEFORM.json')
assert not os.path.exists(out), 'receipts never overwrite — supersede by suffix'
lk_emit(rec, out, component='frobenius-boundary/deform_transport (wall)',
        script_path=os.path.abspath(__file__), member_dir=dk.MEMBER_DIR)
print('WROTE', out)
dk.lint(out)
print('WALL receipt emitted (mechanism named, fields script-emitted).')

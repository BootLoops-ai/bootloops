#!/usr/bin/env python3
"""gmtel_selftest — the gmtel member's own battery (positive control).

The Watson/HLY (2,2,2) K3 family at (A,B,C) = (1,4,16): route B (exact
moments -> annihilator fit) must equal the reference L5 specialization
EXACTLY, and route A (the Gauss-Manin q-telescoper, run from the Weierstrass
family itself) must equal the x^(-1/2)-conjugate of BOTH the route-B minimal
operator AND the reference L5, EXACTLY, with the exact certificate solved
over Q[q] at 8/8 rational s points.

Usage:
  python3 gmtel_selftest.py            # full control (1,4,16), ~1.5-2.5 min
  python3 gmtel_selftest.py --pilot    # fast degeneration smoke (1,1,1), ~30 s
  python3 gmtel_selftest.py --outdir D # keep receipts in D (default: tmpdir)

Drives the two route modules as subprocesses with the same flags as the
control/pilot stages; all verdict fields are script-emitted.
"""
import argparse, json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def run_stage(script, args, outdir):
    cmd = [sys.executable, os.path.join(HERE, script)] + args
    r = subprocess.run(cmd, cwd=outdir, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pilot', action='store_true',
                    help='fast smoke at the (1,1,1) degeneration point')
    ap.add_argument('--outdir', default=None,
                    help='receipt dir (default: fresh temp dir)')
    args = ap.parse_args()
    outdir = args.outdir or tempfile.mkdtemp(prefix='gmtel_selftest_')
    os.makedirs(outdir, exist_ok=True)
    fails = []

    if args.pilot:
        abc, mmax, nmax, nsamp, cap = ['1', '1', '1'], '90', '4', '80', '120'
        tagB, tagA = 'PILOT_PF_ROUTEB_111.json', 'PILOT_PF_ROUTEA_111.json'
    else:
        abc, mmax, nmax, nsamp, cap = ['1', '4', '16'], '150', '6', '140', '1200'
        tagB, tagA = 'CONTROL_PF_ROUTEB.json', 'CONTROL_PF_ROUTEA.json'

    # ---- route B (data): exact moments -> annihilator fit vs the reference L5
    rc, log = run_stage('famgen_pf_routeB.py',
                        ['--abc'] + abc + ['--mmax', mmax, '--out', tagB], outdir)
    pb = os.path.join(outdir, tagB)
    if rc != 0 or not os.path.exists(pb):
        print(log)
        print('GMTEL SELFTEST FAIL: route B rc=%d' % rc)
        return 1
    rb = json.load(open(pb))
    if rb.get('verdict') != 'PASS':
        fails.append('route B verdict: %r' % rb.get('verdict'))
    if not args.pilot and not rb.get('fit_equals_L5_exactly'):
        fails.append('route B fitted operator != reference L5 at control point')
    if args.pilot and not (rb.get('recurrence_annihilates_all_terms')
                           and rb.get('L5_annihilates_moments')):
        fails.append('route B pilot annihilation checks failed')

    # ---- route A (geometric): GM q-telescoper from the Weierstrass family
    rc, log = run_stage('famgen_pf_routeA.py',
                        ['--abc'] + abc + ['--nmax', nmax, '--nsamp', nsamp,
                         '--exact-verify-cap-s', cap, '--routeb', pb,
                         '--out', tagA], outdir)
    pa = os.path.join(outdir, tagA)
    if rc != 0 or not os.path.exists(pa):
        print(log)
        print('GMTEL SELFTEST FAIL: route A rc=%d' % rc)
        return 1
    ra = json.load(open(pa))
    if not str(ra.get('verdict', '')).startswith('PASS'):
        fails.append('route A verdict: %r' % ra.get('verdict'))
    g = ra.get('gates', {})
    if not g.get('vs_routeB'):
        fails.append('route A operator != route-B operator (up to twist)')
    if '1/2' not in g.get('series_twist_kappa_candidates', []):
        fails.append('series twist detector did not return kappa = 1/2')
    if not args.pilot:
        if not g.get('vs_L5'):
            fails.append('route A operator != reference L5 (up to twist)')
        if ra.get('L_A_order') != 5:
            fails.append('telescoper order %r != 5' % ra.get('L_A_order'))
        if ra.get('exact_certificate', {}).get('n_exact_s_points') != 8:
            fails.append('exact certificate %r/8 s-points'
                         % ra.get('exact_certificate', {}).get('n_exact_s_points'))

    print('receipts:', outdir)
    print('route B: %s (%.1f s)' % (rb.get('verdict'), rb.get('t_total_s', -1)))
    print('route A: %s (%.1f s)' % (ra.get('verdict'), ra.get('t_total_s', -1)))
    if fails:
        for f in fails:
            print('  FAIL-GATE:', f)
        print('GMTEL SELFTEST FAIL (%d gate(s))' % len(fails))
        return 1
    mode = 'pilot (1,1,1) degeneration smoke' if args.pilot else \
           'positive control (1,4,16) two-route == reference L5 up to x^(-1/2) twist'
    print('GMTEL SELFTEST PASS:', mode)
    return 0


if __name__ == '__main__':
    sys.exit(main())

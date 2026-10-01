#!/usr/bin/env python3
"""eichler package self-test — ONE command for every member that ships a battery.

    python3 selftest.py [--outdir <scratch>]          (run from tools/eichler/)

Legs, in order:
  M   mirror6 — the positive-control battery (mirror6/run_control.py) run into
      <outdir>/mirror6: pilot then full, 28 boolean rows against published
      instanton numbers (quintic + Rodland, both MUM points); must print
      CONTROL-PASS and exit 0.
  L6  mirror6/run_l6.py — needs the sha-pinned L6_exact.json operator file,
      which is not included: verified to REFUSE by name (fail-closed) when the
      file is absent; reported as a named refusal, never as a silent pass.
  G   genus2 — S1 mestre_port vs the Sage doctest vectors (exact), S2 conic_fast
      worked example with the diagonalization re-checked exactly, R1/R2 the
      fail-closed refusals of invariant_harness / pslq_law_template
      (genus2/selftest.py, imported and run in-process).
  T   theta9 — ships without a battery (its reference receipts are not part of
      this repository): reported as a named SKIP.

Scratch: everything any leg writes goes under --outdir (default: $EICHLER_WORK
if set, else a fresh temporary directory); nothing is written into the package
tree. Exit 0 when M and G pass and L6 refuses as documented; nonzero otherwise.
"""
import argparse
import importlib.util
import os
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def leg_mirror6(work):
    out = os.path.join(work, 'mirror6')
    os.makedirs(out, exist_ok=True)
    t0 = time.monotonic()
    r = subprocess.run([sys.executable, os.path.join(HERE, 'mirror6', 'run_control.py'),
                        '--outdir', out], capture_output=True, text=True, cwd=work)
    msg = (r.stdout or '') + (r.stderr or '')
    tail = '\n'.join(msg.strip().splitlines()[-3:])
    assert r.returncode == 0, f'mirror6 control battery exit {r.returncode}:\n{msg[-800:]}'
    assert 'CONTROL-PASS' in msg, f'mirror6 did not print CONTROL-PASS:\n{msg[-800:]}'
    assert os.path.isfile(os.path.join(out, 'CONTROL_MIRROR6.json')), 'no CONTROL_MIRROR6.json receipt'
    print(tail)
    print(f'M  PASS  mirror6 control battery CONTROL-PASS ({time.monotonic()-t0:.1f}s; receipts in {out})')
    return out


def leg_l6_refusal(ctl_dir):
    missing = os.path.join(ctl_dir, 'L6_exact.json.absent')
    r = subprocess.run([sys.executable, os.path.join(HERE, 'mirror6', 'run_l6.py'),
                        '--l6', missing, '--outdir', ctl_dir],
                       capture_output=True, text=True, cwd=ctl_dir)
    msg = (r.stdout or '') + (r.stderr or '')
    assert r.returncode != 0, 'run_l6 ran without its pinned operator file'
    assert 'L6_exact.json' in msg, f'run_l6 refusal did not name L6_exact.json: {msg[-300:]}'
    print('L6 REFUSED (as documented)  mirror6/run_l6.py names the missing sha-pinned L6_exact.json')


def leg_genus2(work):
    os.environ['EICHLER_WORK'] = os.path.join(work, 'genus2')
    path = os.path.join(HERE, 'genus2', 'selftest.py')
    spec = importlib.util.spec_from_file_location('eichler_genus2_selftest', path)
    g2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g2)          # puts genus2/ on sys.path for its flat imports
    n = g2.run_all()
    print(f'G  PASS  genus2 selftest {n}/{n} legs')


def leg_theta9():
    fx = os.path.join(HERE, 'theta9', 'fixtures')
    nfx = len(os.listdir(fx)) if os.path.isdir(fx) else 0
    print(f'T  SKIP  theta9: no battery ships (reference receipts not distributed; '
          f'{nfx} fixture files present, engine = Eichler.jl via EICHLER_PROJECT)')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n', 1)[0])
    ap.add_argument('--outdir', default=os.environ.get('EICHLER_WORK'),
                    help='scratch/output directory (default: $EICHLER_WORK, else a temp dir)')
    a = ap.parse_args()
    work = os.path.abspath(a.outdir) if a.outdir else tempfile.mkdtemp(prefix='eichler_selftest_')
    os.makedirs(work, exist_ok=True)
    t0 = time.monotonic()
    ctl = leg_mirror6(work)
    leg_l6_refusal(ctl)
    leg_genus2(work)
    leg_theta9()
    print(f'eichler selftest: ALL PASS (M + G; L6 refused by name; T skipped by name) '
          f'wall {time.monotonic()-t0:.1f}s  outdir={work}')


if __name__ == '__main__':
    main()

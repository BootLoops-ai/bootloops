#!/usr/bin/env python3
"""SubTropica battery entry — run: python3 selftest.py

Class: partial (see tools/README.md) — nothing runs end to end without a
HyperFLINT engine build (deps/BUILD.md). What this entry runs:

  R  reference pins: sha256 of every upstream-reference file against
     reference/SHA256SUMS (self-contained; fail loud, names the file).
  F  fibrate control battery: scripts/fibrate/fibrate.py --selfcheck — the
     sound-fibration engine's mandatory controls (PSLQ positive + negative
     controls, adaptive-ZIP-series check vs an exact log truth, GEval/GEvalX
     cross-check). Pure Python (mpmath/sympy), no engine.
  E  driver suite: julia --project=. test/runtests.jl — runs ONLY when
     SUBTROPICA_HF_BIN points at a built engine wrapper (with the GUIDE's
     memory cap); otherwise skips loudly naming the requirement. Building the
     engine and the Julia env is the GUIDE's DEPENDENCIES/ENGINE section.

Exit 0 = R + F pass and E passed or skipped-for-a-named-reason; nonzero
otherwise.
"""
import hashlib
import json
import os
import resource
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MEM_CAP_KB = 32505856          # the GUIDE's ulimit -v for engine and suite


def r_reference_pins():
    sums = os.path.join(HERE, 'reference', 'SHA256SUMS')
    n = 0
    for line in open(sums):
        if not line.strip():
            continue
        want, name = line.split()
        p = os.path.join(HERE, 'reference', name)
        got = hashlib.sha256(open(p, 'rb').read()).hexdigest()
        assert got == want, f'reference/{name}: sha256 {got[:16]}.. != pinned {want[:16]}..'
        n += 1
    print(f'R PASS  reference pins: {n}/{n} files sha256-verified')


def f_fibrate_controls():
    r = subprocess.run([sys.executable,
                        os.path.join(HERE, 'scripts', 'fibrate', 'fibrate.py'),
                        '--selfcheck'], capture_output=True, text=True)
    assert r.returncode == 0, f'fibrate --selfcheck rc={r.returncode}: {r.stderr[-400:]}'
    rep = json.loads(r.stdout)
    assert rep.get('selfcheck') == 'PASS', rep
    print(f'F PASS  fibrate control battery: {rep}')


def e_driver_suite():
    hf = os.environ.get('SUBTROPICA_HF_BIN')
    if not hf or not os.path.isfile(hf):
        print('E SKIP  driver suite: SUBTROPICA_HF_BIN is not set (or not a '
              'file) — the Julia suite needs a built HyperFLINT engine '
              '(deps/BUILD.md) plus SUBTROPICA_MZV_DATA and the instantiated '
              'Julia project; see GUIDE.md. Skipping is the designed '
              'engine-absent behavior, not breakage.')
        return 0
    def cap():
        resource.setrlimit(resource.RLIMIT_AS,
                           (MEM_CAP_KB * 1024, MEM_CAP_KB * 1024))
    r = subprocess.run(['julia', '--project=' + HERE,
                        os.path.join(HERE, 'test', 'runtests.jl')],
                       cwd=HERE, preexec_fn=cap)
    if r.returncode == 0:
        print('E PASS  driver suite (test/runtests.jl) green')
    else:
        print(f'E FAIL  driver suite rc={r.returncode}', file=sys.stderr)
    return r.returncode


if __name__ == '__main__':
    r_reference_pins()
    f_fibrate_controls()
    rc = e_driver_suite()
    if rc == 0:
        print('subtropica selftest: PASS (R, F; E per engine availability above)')
    sys.exit(rc)

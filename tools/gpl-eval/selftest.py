#!/usr/bin/env python3
"""gpl-eval battery entry — run: python3 selftest.py

Thin wrapper over the package's own test suite (GPLEval.jl/test/runtests.jl),
honoring the partial-class contract: a missing engine or dependency is a loud
NAMED skip, never a raw load traceback.

  1. no `julia` on PATH -> named SKIP (exit 0). GPLEval.jl is a Julia package
     (julia >= 1.10); install julia (e.g. via juliaup), then rerun.
  2. cheap environment probe (package lookup, nothing loaded). If the project
     cannot resolve its dependencies yet, wire it up: first
     `Pkg.instantiate()` against the committed Project/Manifest, then the full
     GPLEval.jl/deps/setup.jl wiring (develops the in-repo Eichler +
     LandauAlphabet packages, instantiates the registered deps Arblib + JSON3
     from the General registry — needs network on a cold machine). If the
     wiring cannot complete -> named SKIP (exit 0) stating the exact commands
     and package names.
  3. run the suite: julia --project=GPLEval.jl GPLEval.jl/test/runtests.jl.
     The GiNaC cross-check legs skip with a warning when the external GiNaC
     C++ library is unavailable; that skip is by design (class: partial), not
     breakage.

Exit 0 = suite green, or a named dependency skip; nonzero = a genuine test
failure, or a fail-closed refusal (named) when the environment resolves but
fails to load.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, 'GPLEval.jl')

# Package lookup only — nothing is loaded, so this is fast and cannot traceback.
PROBE = ('for name in ("GPLEval", "Arblib", "JSON3", "Eichler", "LandauAlphabet")\n'
         '    id = Base.identify_package(name)\n'
         '    if id === nothing || Base.locate_package(id) === nothing\n'
         '        print(name); exit(3)\n'
         '    end\n'
         'end\n')

LOAD_FAIL_SIGNS = ('not found in current path',
                   'does not seem to be installed',
                   'Failed to precompile',
                   'error while loading shared libraries')


def probe(julia):
    r = subprocess.run([julia, '--startup-file=no', f'--project={PKG}',
                        '-e', PROBE], capture_output=True, text=True)
    return r.returncode == 0, r.stdout.strip()


def named_skip(reason):
    print(f'''
gpl-eval selftest: SKIPPED (named dependency gate) — {reason}
GPLEval.jl needs the Julia packages Arblib + JSON3 (General registry) and the
two packages that ship in this repository, Eichler (upgrades/Eichler.jl) and
LandauAlphabet (tools/landau-alphabet/LandauAlphabet.jl). To wire a cold
clone (one time, needs network for the registered deps):

    julia --project=GPLEval.jl -e 'import Pkg; Pkg.instantiate()'
    # or, if that fails (e.g. a different julia version than the committed
    # Manifest), the full re-wiring:
    julia GPLEval.jl/deps/setup.jl

then rerun:  python3 selftest.py
This skip is the designed dependency-absent behavior for a partial-class
package, not breakage. Nothing was verified.''')
    return 0


def main():
    julia = shutil.which('julia')
    if not julia:
        print('gpl-eval selftest: SKIPPED (named engine gate) — no `julia` on '
              'PATH. GPLEval.jl is a Julia package (julia >= 1.10); install '
              'julia (e.g. via juliaup), then rerun. This skip is the designed '
              'engine-absent behavior for a partial-class package, not '
              'breakage. Nothing was verified.')
        return 0

    ok, missing = probe(julia)
    if not ok:
        print(f'gpl-eval selftest: environment not wired yet (first missing: '
              f'{missing or "unknown"}) — wiring now (downloads the registered '
              f'deps Arblib + JSON3 on a cold machine; may take a few minutes) '
              f'...', flush=True)
        # fast path: the committed Project/Manifest (in-repo dev paths are
        # relative, so a fresh clone resolves them as-is)
        subprocess.run([julia, '--startup-file=no', f'--project={PKG}',
                        '-e', 'import Pkg; Pkg.instantiate()'], cwd=HERE)
        ok, missing = probe(julia)
    if not ok:
        # full re-wiring (re-resolves for this julia version)
        r = subprocess.run([julia, os.path.join(PKG, 'deps', 'setup.jl')],
                           cwd=HERE)
        if r.returncode == 0:
            ok, missing = probe(julia)
    if not ok:
        return named_skip(f'the Julia environment could not be wired '
                          f'(first unresolved package: {missing or "unknown"})')

    r = subprocess.run([julia, f'--project={PKG}',
                        os.path.join(PKG, 'test', 'runtests.jl')],
                       cwd=HERE, capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    sys.stderr.write(r.stderr)
    if r.returncode != 0:
        blob = r.stdout + r.stderr
        if any(sig in blob for sig in LOAD_FAIL_SIGNS):
            print('gpl-eval selftest: REFUSED (fail-closed) — the Julia '
                  'environment resolved but failed to LOAD (this is not a '
                  'test failure). Re-wire it with `julia '
                  'GPLEval.jl/deps/setup.jl` and rerun.', file=sys.stderr)
            return 1
    return r.returncode


if __name__ == '__main__':
    sys.exit(main())

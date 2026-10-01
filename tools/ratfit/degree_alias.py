#!/usr/bin/env python3
"""ratfit.degree_alias — fit-degree probe: grid-alias vs height discriminator
(member of the ratfit package; law receipt DEGREE_ALIAS_DISCRIMINATOR.json
ships as package data beside this file, two-sided s5806/s5805 oracle proof).

SYMPTOM: CRT +prime passes NEVER lift rational-reconstruction residuals
(MQRR failures persist pass after pass). A position-only classifier (interior
vs grid-top sample position) is BLIND to grid aliasing — an interpolation
alias is consistent across ALL primes and MQRR reads it as interior, so it
mislabels as HEIGHT and burns +prime passes with 0 lift (s5806: 0/272 lifted
at 8p AND 9p).

THE LAW (per-prime fit degrees on the FAILING entries decide the failure axis):
- num degrees PINNED AT THE GRID CEILING (nd-1 / neta-1) with TRIVIAL
  denominator (0,0), identical across all probed primes
    => DEGREE-ALIASED (grid too small) -> fix: single-axis grid EXTENSION
       (never +prime).  Proof side A: s5806, grid (50,100), (49,99)/(0,0)
       identical at 7 primes; corroboration 0/272 lifted at 8p and 9p.
- degrees stable BELOW ceiling with NONTRIVIAL denominator
    => HEIGHT (CRT modulus too small) -> fix: +prime passes (never grid).
       Proof side B: s5805, grid (29,101), num (11-13,45-47) den (9-10,46-48)
       at 5 primes; corroboration: all 1031 lifted by +1 prime.
- cross-prime degree DISAGREEMENT => INCONCLUSIVE: probe more primes/cells;
  never pick a fix on this evidence (fail-closed).

IMPORT: `from ratfit.degree_alias import classify, classify_cell`
(or `import degree_alias` with tools/ratfit on sys.path — the module has no
dependencies beyond the standard library).

API (pure functions, no I/O):
  classify_cell(grid, fits) -> 'ALIAS' | 'HEIGHT' | 'INCONCLUSIVE'
      grid = (nd, ne) point counts per axis;
      fits = per-prime list of ((num_d, num_e), (den_d, den_e)) fit degrees
             for ONE failing cell/coefficient.
  classify(grid, cells) -> {'verdict', 'axis', 'counts'}
      cells = list of per-cell fits lists. verdict ALIAS only when EVERY
      fitted cell is alias-pinned (the reference n_alias >= n_fits rule);
      axis = 'd' if #d-pinned >= #e-pinned (reference tie-break), else 'eta'.
  discriminator() -> dict
      the law receipt (DEGREE_ALIAS_DISCRIMINATOR.json, package data) parsed;
      DISCRIMINATOR_PATH is its absolute path.
CLI: python3 -m ratfit.degree_alias --selftest      (from a dir with tools/ on
     PYTHONPATH), or python3 tools/ratfit/degree_alias.py --selftest
     (replays the two-sided proof + two negatives + the receipt cross-check).
"""
import json
import os

DISCRIMINATOR_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  'DEGREE_ALIAS_DISCRIMINATOR.json')


def discriminator(path=None):
    """Return the parsed law receipt shipped beside this module."""
    with open(path or DISCRIMINATOR_PATH) as fh:
        return json.load(fh)


def classify_cell(grid, fits):
    nd, ne = grid
    if not fits:
        return 'INCONCLUSIVE'
    if any(f != fits[0] for f in fits[1:]):
        return 'INCONCLUSIVE'          # cross-prime disagreement: probe more
    (num_d, num_e), (den_d, den_e) = fits[0]
    den_trivial = (den_d, den_e) == (0, 0)
    pinned_d = num_d >= nd - 1
    pinned_e = num_e >= ne - 1
    if den_trivial and (pinned_d or pinned_e):
        return 'ALIAS'
    if not den_trivial and num_d < nd - 1 and num_e < ne - 1 \
            and den_d < nd - 1 and den_e < ne - 1:
        return 'HEIGHT'
    return 'INCONCLUSIVE'


def classify(grid, cells):
    nd, ne = grid
    n_alias = n_height = n_inc = n_alias_d = n_alias_e = 0
    for fits in cells:
        c = classify_cell(grid, fits)
        if c == 'ALIAS':
            n_alias += 1
            (num_d, num_e), _ = fits[0]
            if num_d >= nd - 1:
                n_alias_d += 1
            if num_e >= ne - 1:
                n_alias_e += 1
        elif c == 'HEIGHT':
            n_height += 1
        else:
            n_inc += 1
    n_fits = n_alias + n_height + n_inc
    counts = {'n_fits': n_fits, 'alias': n_alias, 'height': n_height,
              'inconclusive': n_inc, 'alias_d': n_alias_d, 'alias_e': n_alias_e}
    if n_fits == 0 or n_inc > 0:
        return {'verdict': 'INCONCLUSIVE', 'axis': None, 'counts': counts}
    if n_alias >= n_fits:               # reference rule: ALL fitted cells pinned
        axis = 'd' if n_alias_d >= n_alias_e else 'eta'
        return {'verdict': 'ALIAS', 'axis': axis, 'counts': counts}
    return {'verdict': 'HEIGHT', 'axis': None, 'counts': counts}


def _selftest():
    # side A — s5806 shape: grid (50,100), (49,99)/(0,0) identical at 7 primes
    a = classify((50, 100), [[((49, 99), (0, 0))] * 7] * 5)
    assert a['verdict'] == 'ALIAS' and a['axis'] == 'd', a
    # side B — s5805 shape: grid (29,101), below-ceiling num + nontrivial den, 5 primes
    b = classify((29, 101), [[((12, 46), (9, 47))] * 5,
                             [((13, 45), (10, 48))] * 5,
                             [((11, 47), (9, 46))] * 5])
    assert b['verdict'] == 'HEIGHT', b
    # negative — cross-prime disagreement must refuse, never pick a fix
    c = classify((50, 100), [[((49, 99), (0, 0)), ((30, 99), (0, 0))]])
    assert c['verdict'] == 'INCONCLUSIVE', c
    # negative — pinned num but NONTRIVIAL den is not the alias signature
    d = classify((50, 100), [[((49, 20), (3, 2))] * 5])
    assert d['verdict'] == 'INCONCLUSIVE', d
    # receipt cross-check — the shipped law receipt must carry the two proof
    # sides replayed above (grids of side A / side B), so code and receipt
    # cannot drift apart silently.
    r = discriminator()
    assert tuple(r['side_A_alias']['grid']) == (50, 100), r['side_A_alias']
    assert tuple(r['side_B_height']['grid']) == (29, 101), r['side_B_height']
    assert 'law' in r and r['law'], r
    print('degree_alias selftest PASS (A=ALIAS/d, B=HEIGHT, negatives refuse, '
          'receipt grids match)')


if __name__ == '__main__':
    import sys
    if '--selftest' in sys.argv:
        _selftest()
    else:
        print(__doc__)

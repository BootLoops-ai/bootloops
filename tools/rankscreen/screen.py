#!/usr/bin/env python3
"""Multi-prime parallel rank screen — the two-precision law applied to
linear algebra.

Runs the mod-p RREF (modp_rref.py, verdict semantics of the certified exact
eliminators) for k primes IN PARALLEL (one process per prime — farm-unit
law), then compares (rank, pivot sequence, inconsistent-row set) across
primes:

  identical across all k primes  -> CERTIFIED-SCREEN
        + closure_candidate flag (rank-full on the watched column block AND
          zero inconsistent rows) -> only then run the exact-Q solve,
          seeded with the agreed pivot structure (shim.exact_seeded).
  any disagreement               -> ESCALATE-TO-EXACT (never guess).

A CERTIFIED-SCREEN closure candidate is NEVER promoted to a closure by the
screen alone: the exact solve is always the closure authority.  What the
screen certifies cheaply is the REFUTATION side (rank deficits and
inconsistency sets) and the pivot structure.  False-agreement probability
bounds: see GUIDE.md, stated honestly —
any wrong unanimous verdict forces prod(p_i) to divide a specific nonzero
minor-class integer invariant of the augmented system.

CLI:
  python3 screen.py ROWS.json.gz [--k 2] [--backend auto|dense|sparse]
                    [--nunk N] [--out RECEIPT.json] [--primes p1,p2,...]
ROWS.json.gz = rankscreen-rows-v1 (shim.serialize_rows).

API: screen_rows(rows, labels, ncols, ...) -> verdict dict (same content as
the receipt).  Uses fork + module globals: the parent's row list is shared
copy-on-write with the per-prime workers — rows are never pickled.
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modp_rref import (BadPrime, PRIMES_DENSE, PRIMES_SPARSE, rref_modp,
                       DENSE_NCOLS_MAX)

_G = {}  # fork-shared: rows, ncols, backend


def _worker(p):
    try:
        r = rref_modp(_G['rows'], p, _G['ncols'], _G['backend'])
        return {'prime': p, 'ok': True, 'rank': r.rank,
                'pivots': list(r.pivots), 'incon_idx': list(r.incon_idx),
                'backend': r.backend, 'secs': r.secs}
    except BadPrime as e:
        return {'prime': p, 'ok': False, 'badprime': str(e)}


def _pool(backend, ncols):
    if backend == 'sparse':
        return PRIMES_SPARSE
    if backend == 'dense':
        return PRIMES_DENSE
    return (PRIMES_DENSE if (ncols is not None and ncols <= DENSE_NCOLS_MAX)
            else PRIMES_SPARSE)


def screen_rows(rows, labels=None, ncols=None, k=2, backend='auto',
                primes=None, nunk=None, receipt_path=None, rows_sha=None,
                nproc=None, tag=''):
    """rows: [({col: Fraction}, Fraction rhs)]; labels parallel list or None.
    Returns the verdict dict (and writes it to receipt_path if given)."""
    assert k >= 2, "screen requires >= 2 primes (two-precision law)"
    t0 = time.time()
    if ncols is None:
        ncols = 1 + max((max(d) for d, _ in rows if d), default=-1)
    pool = list(primes) if primes else list(_pool(backend, ncols))
    assert len(pool) >= k, "prime pool exhausted"
    _G['rows'] = rows; _G['ncols'] = ncols; _G['backend'] = backend

    ctx = mp.get_context('fork')
    results, replaced, queue = [], [], pool[:]
    active = queue[:k]; queue = queue[k:]
    while True:
        with ctx.Pool(processes=min(len(active), nproc or len(active))) as pl:
            got = pl.map(_worker, active)
        bad = [g for g in got if not g['ok']]
        results += [g for g in got if g['ok']]
        if not bad:
            break
        for b in bad:
            replaced.append({'prime': b['prime'], 'why': b['badprime']})
        if len(queue) < len(bad):
            raise BadPrime(f"prime pool exhausted; replaced={replaced}")
        active = queue[:len(bad)]; queue = queue[len(bad):]
    results.sort(key=lambda g: g['prime'], reverse=True)

    sig0 = (results[0]['rank'], results[0]['pivots'], results[0]['incon_idx'])
    agree = all((g['rank'], g['pivots'], g['incon_idx']) == sig0
                for g in results[1:])
    verdict = 'CERTIFIED-SCREEN' if agree else 'ESCALATE-TO-EXACT'
    out = {
        'tool': 'rankscreen', 'format': 'rankscreen-verdict-v1', 'tag': tag,
        'ncols': ncols, 'n_rows': len(rows), 'k': k,
        'primes': [g['prime'] for g in results],
        'primes_replaced': replaced,
        'backend': results[0]['backend'],
        'verdict': verdict, 'agree': agree,
        'rank': sig0[0] if agree else None,
        'n_inconsistent': len(sig0[2]) if agree else None,
        'per_prime': [{'prime': g['prime'], 'rank': g['rank'],
                       'n_inconsistent': len(g['incon_idx']),
                       'pivots_sha': hashlib.sha256(
                           json.dumps(g['pivots']).encode()).hexdigest()[:16],
                       'secs': g['secs']} for g in results],
        'wall_secs': round(time.time() - t0, 3),
    }
    if rows_sha:
        out['rows_sha256'] = rows_sha
    if agree:
        piv = sig0[1]
        out['pivot_cols'] = [c for c, _ in piv]
        out['pivot_rows'] = [ri for _, ri in piv]
        out['incon_idx'] = list(sig0[2])
        if labels is not None:
            out['incon_labels'] = [labels[i] for i in sig0[2]]
        if nunk is not None:
            npure = sum(1 for c, _ in piv if c < nunk)
            out['nunk'] = nunk
            out['pure_pivots'] = npure
            out['closure_candidate'] = (npure == nunk
                                        and len(sig0[2]) == 0)
        else:
            out['closure_candidate'] = (sig0[0] == ncols
                                        and len(sig0[2]) == 0)
    else:
        out['disagreement'] = [
            {'prime': g['prime'], 'rank': g['rank'],
             'n_inconsistent': len(g['incon_idx'])} for g in results]
    if receipt_path:
        json.dump(out, open(receipt_path, 'w'), indent=1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rows_file')
    ap.add_argument('--k', type=int, default=2)
    ap.add_argument('--backend', default='auto',
                    choices=['auto', 'dense', 'sparse'])
    ap.add_argument('--nunk', type=int, default=None,
                    help='pure-column count for the closure-candidate flag')
    ap.add_argument('--primes', default=None,
                    help='comma-separated prime override')
    ap.add_argument('--out', default=None)
    ap.add_argument('--tag', default='')
    a = ap.parse_args()
    from shim import load_rows, sha256_file
    rows, labels, ncols, meta = load_rows(a.rows_file)
    primes = ([int(x) for x in a.primes.split(',')] if a.primes else None)
    out = screen_rows(rows, labels, ncols, k=a.k, backend=a.backend,
                      primes=primes, nunk=a.nunk, receipt_path=a.out,
                      rows_sha=sha256_file(a.rows_file), tag=a.tag)
    print(json.dumps({k2: v for k2, v in out.items()
                      if k2 not in ('pivot_cols', 'pivot_rows', 'incon_idx',
                                    'incon_labels')}, indent=1))
    print(f"VERDICT: {out['verdict']}"
          + (f" rank {out['rank']}/{ncols} inc {out['n_inconsistent']}"
             f" closure_candidate={out.get('closure_candidate')}"
             if out['agree'] else "  (primes disagree — run exact)"))


if __name__ == '__main__':
    main()

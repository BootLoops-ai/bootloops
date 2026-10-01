#!/usr/bin/env python3
"""Integration shim for the pm_core / pm_ext closure-engine lineage.

The closure engines of the gauss_rank_mpq class build
    rows_all = [(label, {col: Fraction}, Fraction rhs)]
and call an exact eliminator (elim.gauss_rank_mpq / window_core.gauss_rank)
returning (rank, ninc, piv, red, redr, incon_labels).  This shim gives them:

  serialize_rows / load_rows   the standard on-disk form (rankscreen-rows-v1)
  screen_rows_labeled          screen straight off the in-memory rows_all
  gauss_rank_screened          DROP-IN for gauss_rank_mpq at the full-system
                               call site: screens first; runs the exact-Q
                               solve ONLY on closure candidates (or on prime
                               disagreement), seeded with the mod-p pivot
                               structure; returns the same 6-tuple shape.
  exact_seeded                 exact mpq elimination with the screen's pivot
                               plan (pivot rows first, prediction verified,
                               automatic fallback to plain order).

DO NOT retrofit into a live running chain — adopt at the next natural run
(edit the elimination call site; everything upstream is untouched).

Return-shape note for the screened (non-exact) branch:
  red, redr are None (no exact reduced rows exist);
  piv is {pivot_col: original_row_index} — key-compatible with the engines'
  pivot-count/membership uses (`[c for c in piv if c < NUNK]` etc.).
The closure branch (and any disagreement) always returns the genuine exact
6-tuple, so downstream solution extraction is unchanged.
"""
import gzip
import hashlib
import json
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from screen import screen_rows

FORMAT = 'rankscreen-rows-v1'


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for blk in iter(lambda: f.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest()


# ------------------------------------------------------------- serialize ----
def serialize_rows(rows_all, ncols, path, meta=None):
    """rows_all: [(label, {col: Fraction}, Fraction rhs)] -> gzip JSON."""
    doc = {'format': FORMAT, 'ncols': ncols, 'meta': meta or {},
           'labels': [str(lab) for lab, _, _ in rows_all],
           'rows': [[{str(c): str(q) for c, q in cv.items()}, str(rhs)]
                    for _, cv, rhs in rows_all]}
    with gzip.open(path, 'wt') as f:
        json.dump(doc, f)
    return sha256_file(path)


def load_rows(path):
    """-> (rows [( {col:Fraction}, Fraction )], labels, ncols, meta)."""
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        doc = json.load(f)
    assert doc['format'] == FORMAT, doc.get('format')
    rows = [({int(c): Fraction(q) for c, q in cv.items()}, Fraction(rhs))
            for cv, rhs in doc['rows']]
    return rows, doc['labels'], doc['ncols'], doc.get('meta', {})


# ------------------------------------------------------- exact reference ----
def gauss_rank_mpq(rowlist, track_labels=None, progress=0, tag=''):
    """Exact-Q eliminator, gmpy2.mpq — ALGORITHM-IDENTICAL to the engines'
    certified elim.gauss_rank_mpq / window_core.gauss_rank (same row order,
    min-column pivot, insertion-order reduction).  Kept here so the tool has
    no import edge into a run directory; gate-verified to reproduce the
    reference receipts on the real systems."""
    from gmpy2 import mpq
    piv = {}; red = []; redr = []; rank = 0; incon = []
    labels = track_labels if track_labels is not None \
        else [None] * len(rowlist)
    t0 = time.time()
    for idx, ((rdict0, rv0), lab) in enumerate(zip(rowlist, labels)):
        rdict = {c: mpq(v) for c, v in rdict0.items() if v != 0}
        rv = mpq(rv0)
        for col, ri in list(piv.items()):
            if col in rdict and rdict[col] != 0:
                f = rdict[col] / red[ri][col]
                for c2, v2 in red[ri].items():
                    rdict[c2] = rdict.get(c2, 0) - f * v2
                    if rdict[c2] == 0: del rdict[c2]
                rv -= f * redr[ri]
        rdict = {k: v for k, v in rdict.items() if v != 0}
        if rdict:
            c0 = min(rdict)
            piv[c0] = len(red); red.append(rdict); redr.append(rv)
            rank += 1
        elif rv != 0:
            incon.append(lab)
        if progress and idx and idx % progress == 0:
            print(f"  elim{tag}: row {idx}/{len(rowlist)} rank {rank} "
                  f"inc {len(incon)} ({time.time()-t0:.0f}s)", flush=True)
    return rank, len(incon), piv, red, redr, incon


# ---------------------------------------------------------- seeded exact ----
def exact_seeded(rowlist, track_labels, pivot_rows, progress=0, tag=''):
    """Exact solve seeded with the screen's agreed pivot structure:
    process the predicted pivot rows FIRST (original relative order — this
    reproduces the same pivot basis when the prediction is exact-correct),
    then reduce every remaining row against the finished basis (pure
    verification sweep).  Rank / pivot SET / inconsistent-row SET are
    order-invariant given the same pivot-row set, so results are the
    engines' own (property + gate-verified).  If the prediction breaks
    (a predicted pivot row fails to yield a pivot, or a predicted free row
    yields one), falls back to the plain-order exact eliminator."""
    labels = track_labels if track_labels is not None \
        else [None] * len(rowlist)
    pset = set(pivot_rows)
    order = [i for i in pivot_rows] + \
            [i for i in range(len(rowlist)) if i not in pset]
    rl = [rowlist[i] for i in order]
    lb = [labels[i] for i in order]
    rank, ninc, piv, red, redr, incon = gauss_rank_mpq(
        rl, lb, progress=progress, tag=tag + ':seeded')
    npred = len(pivot_rows)
    piv_from_pred = sum(1 for c, ri in piv.items() if ri < npred)
    if rank != npred or piv_from_pred != npred:
        print(f"  rankscreen: seed prediction broke (rank {rank} vs "
              f"predicted {npred}) — falling back to plain order",
              flush=True)
        return gauss_rank_mpq(rowlist, track_labels, progress=progress,
                              tag=tag + ':fallback')
    return rank, ninc, piv, red, redr, incon


# --------------------------------------------------------------- drop-in ----
def screen_rows_labeled(rows_all, ncols, nunk=None, k=2, backend='auto',
                        receipt_path=None, tag=''):
    """rows_all = [(label, cv, rhs)] straight from the engine."""
    labels = [lab for lab, _, _ in rows_all]
    rows = [(cv, rhs) for _, cv, rhs in rows_all]
    return screen_rows(rows, labels, ncols, k=k, backend=backend,
                       nunk=nunk, receipt_path=receipt_path, tag=tag)


def gauss_rank_screened(rowlist, track_labels=None, progress=0, tag='',
                        ncols=None, nunk=None, k=2, backend='auto',
                        receipt_path=None, force_exact=False):
    """Drop-in for elim.gauss_rank_mpq(rowlist, labels, progress, tag).

    Screens with k primes in parallel; runs the exact-Q solve ONLY when the
    screen certifies a closure candidate (then seeded with the agreed pivot
    structure) or when the primes disagree.  Otherwise returns the
    mod-p-certified refutation verdict in the engine's 6-tuple shape with
    red = redr = None (see module docstring)."""
    v = screen_rows(rowlist, track_labels, ncols, k=k, backend=backend,
                    nunk=nunk, receipt_path=receipt_path, tag=tag)
    if not v['agree']:
        print(f"  rankscreen[{tag}]: ESCALATE-TO-EXACT "
              f"({v['disagreement']})", flush=True)
        return gauss_rank_mpq(rowlist, track_labels, progress, tag)
    print(f"  rankscreen[{tag}]: CERTIFIED-SCREEN rank {v['rank']}/"
          f"{v['ncols']} inc {v['n_inconsistent']} "
          f"closure_candidate={v.get('closure_candidate')} "
          f"({v['wall_secs']}s, primes {v['primes']})", flush=True)
    if force_exact or v.get('closure_candidate'):
        return exact_seeded(rowlist, track_labels, v['pivot_rows'],
                            progress=progress, tag=tag)
    labels = track_labels if track_labels is not None \
        else [None] * len(rowlist)
    piv = {c: ri for c, ri in zip(v['pivot_cols'], v['pivot_rows'])}
    incon = [labels[i] for i in v['incon_idx']]
    return v['rank'], v['n_inconsistent'], piv, None, None, incon

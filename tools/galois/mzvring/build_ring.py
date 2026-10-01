#!/usr/bin/env python3
# galois member — relation-ring builder CLI.
"""build_ring.py — build the exact MZV relation ring, pickle the tables.

Path-robust — no hardcoded paths; explicit --wmax and
--out. Finished builds (ring17/ring21/ring23.pkl) go to whatever bank
directory you keep (see package __init__ RING_BANK).

Usage: python3 tools/galois/mzvring/build_ring.py [--wmax N] [--out /path/ringN.pkl]
                             [--resume PKL] [--checkpoint]
Default --out: ./ring<N>.pkl in the current directory.
--resume PKL: load an existing table pkl and continue from its wmax+1 (the
  w-elimination at every NEW weight still runs in full; only lower-weight
  canonicalization tables are reused — they are exactly what a fresh build
  would have produced, and the bank pkls carry their own validation receipts).
--checkpoint: atomically dump <out>.partial after EVERY weight — an
  externally killed build otherwise loses hours of state to a clean log.
Cost (measured, single core, one recent x86 host): wmax=21 ~30 min, wmax=23 ~2 h.
"""
import argparse
import os
import pickle
import resource
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
from galois.mzvring import mzvring as MR  # the member package, by identity


def build(wmax, out_path, resume=None, checkpoint=False):
    t0 = time.time()
    R = MR.Ring.__new__(MR.Ring)
    R.canon, R.survivors, R.wmax = {}, {}, wmax
    w0 = 3
    if resume:
        with open(resume, 'rb') as f:
            T = pickle.load(f)
        R.canon, R.survivors = T['canon'], T['survivors']
        w0 = T['wmax'] + 1
        print(f"RESUME from {resume} (wmax {T['wmax']}) -> building "
              f"{w0}..{wmax}", flush=True)
    for w in range(w0, wmax + 1):
        t1 = time.time()
        R._build(w)
        print(f"w={w}: {time.time()-t1:6.1f}s  survivors {R.survivors[w]}",
              flush=True)
        if checkpoint:
            tmp = out_path + '.partial.tmp'
            with open(tmp, 'wb') as f:
                pickle.dump({'canon': R.canon, 'survivors': R.survivors,
                             'wmax': w}, f)
            os.replace(tmp, out_path + '.partial')
    with open(out_path, 'wb') as f:
        pickle.dump({'canon': R.canon, 'survivors': R.survivors,
                     'wmax': wmax}, f)
    rss_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == 'darwin':
        rss_kib //= 1024  # ru_maxrss is bytes on macOS, KiB on Linux
    print(f"DONE {time.time()-t0:.0f}s RSS_MB {rss_kib//1024} "
          f"-> {out_path}", flush=True)
    return R


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--wmax', type=int, default=21)
    ap.add_argument('--out', default=None)
    ap.add_argument('--resume', default=None)
    ap.add_argument('--checkpoint', action='store_true')
    a = ap.parse_args()
    build(a.wmax, a.out or f'ring{a.wmax}.pkl', a.resume, a.checkpoint)

#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""memory_basis.py -- the memory basis of a longhand.disteval run BY OBJECT: the figure a caller declares on the fence (memory.max of the cgroup
leaf) for W disteval workers with a resident python driver.
  formula (the nominal per-worker form):  memory_max_bytes(W) = page_up(1.3 x W x per_worker_maxrss_B)
  corrected (the form USED):            memory_max_bytes(W) = page_up(1.3 x (driver_maxrss_B + W x per_worker_maxrss_B))
The corrected form is the one used, because the disteval workers are ~6 MB C++ processes and the python driver carries ~90 MB, so the
per-worker formula alone leaves no room for the driver and would OOM at spawn (driver_fits_in_formula_headroom False in the record).  The W 8
PROBE label doubles the basis once (memory.max = 2 x basis); the raise-once knob is 2 x memory.max.  The measured figures of record travel in
../examples/ndpent_top/BASIS_record.json (per-worker maxrss = max VmHWM over the pysecdec_cpuworker pids of a two-worker smoke sampled every 0.2 s; the
driver's VmHWM beside; GNU-time maxrss and the leaf memory.peak beside).  page_up and maxrss_kb are the chain's units of record, unchanged.
usage: memory_basis.py --W 8 [--basis-file BASIS.json | --worker-bytes N --driver-bytes N] [--probe-W 8]"""
import os, sys, json, re, subprocess, hashlib, math, argparse


def utc():
    return subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'], capture_output=True, text=True).stdout.strip()


def sha256(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def page_up(b):
    return ((int(b) + 4095) // 4096) * 4096


def maxrss_kb(timefile):
    m = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', open(timefile).read())
    return int(m.group(1)) if m else None


def basis_entry(W, worker_B, driver_B, probe_W=8):
    """one per-W entry in the record's form: the nominal formula, the headroom test, the corrected figure, the PROBE doubling, the raise-once."""
    formula = page_up(1.3 * W * worker_B)
    corrected = page_up(1.3 * (driver_B + W * worker_B)) if driver_B else None
    probe = (W == probe_W)
    base = corrected if corrected else formula
    cap = 2 * base if probe else base
    headroom = formula - W * worker_B
    return {'W': W, 'per_worker_maxrss_B': worker_B, 'W_x_worker_B': W * worker_B, 'driver_maxrss_B': driver_B,
            'formula_memory_max_B(per_worker x W x 1.3, page_up)': formula, 'formula_headroom_over_W_x_worker_B': headroom, 'driver_fits_in_formula_headroom': (driver_B is not None and driver_B <= headroom),
            'corrected_memory_max_B(1.3 x (driver + W x per_worker), page_up)': corrected,
            'label': 'PROBE' if probe else 'UNCONTENDED', 'memory_max_bytes': cap, 'memory_max_GiB': round(cap / 2 ** 30, 3),
            'memory_max_rule': ('W %d = PROBE: memory.max = the labeled 2x of the W-%d basis' % (W, W) if probe else 'memory.max = the basis') + '; the basis USED = the CORRECTED figure (driver + W x per-worker) x 1.3 -- the nominal per-worker formula alone undershoots (see driver_fits_in_formula_headroom)',
            'raise_once_bytes': 2 * cap, 'raise_once_rule': 'ONE labeled 2x raise at 0.85 x memory.max, rewriting the declaration file', 'census_memg_gib': math.ceil(cap / 2 ** 30)}


def read_basis_file(path):
    """the measured figures from a BASIS json (source.per_worker_maxss.bytes / source.driver_maxrss.bytes) or a smoke receipt (per_worker_maxrss / driver_maxrss)."""
    R = json.load(open(path))
    src = R.get('source', R)
    worker_B = int(src['per_worker_maxrss']['bytes'])
    driver_B = int(src['driver_maxrss']['bytes']) if src.get('driver_maxrss', {}).get('bytes') else None
    return worker_B, driver_B, R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--W', type=int, required=True)
    ap.add_argument('--basis-file', default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'examples', 'ndpent_top', 'BASIS_record.json'))
    ap.add_argument('--worker-bytes', type=int, default=None)
    ap.add_argument('--driver-bytes', type=int, default=None)
    ap.add_argument('--probe-W', type=int, default=8)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.worker_bytes:
        worker_B, driver_B, src = a.worker_bytes, a.driver_bytes, {'given': 'command line'}
        kind = 'GIVEN on the command line (declare its measurement beside)'
    else:
        if not os.path.exists(a.basis_file):
            print('REFUSED: basis file %s absent -- the basis is MEASURED or given (--worker-bytes/--driver-bytes); no declared basis is emitted' % a.basis_file)
            return 2
        worker_B, driver_B, src = read_basis_file(a.basis_file)
        kind = 'MEASURED (%s sha16 %s)' % (os.path.basename(a.basis_file), sha256(a.basis_file)[:16])
    e = basis_entry(a.W, worker_B, driver_B, a.probe_W)
    out = {'receipt': 'MEMORY BASIS by object: memory.max(W) = page_up(1.3 x (driver + W x per_worker)); PROBE W = 2x; raise-once = 2x memory.max', 'stamp_utc': utc(), 'basis_kind': kind, 'entry': e,
           'declare': 'memory.max = %d B (%.3f GiB) for W %d [%s]; cpu.max = %d 100000 (= W x 100000); swap 0; the raise-once knob %d B' % (e['memory_max_bytes'], e['memory_max_GiB'], a.W, e['label'], a.W * 100000, e['raise_once_bytes'])}
    if a.out:
        json.dump(out, open(a.out, 'w'), indent=1)
    print('BASIS W=%d %s: memory.max %d B (%.3f GiB) [nominal formula %d B, headroom %d B, driver %s B fits=%s; corrected %s B]; raise-once %d B; basis %s' % (a.W, e['label'], e['memory_max_bytes'], e['memory_max_GiB'], e['formula_memory_max_B(per_worker x W x 1.3, page_up)'], e['formula_headroom_over_W_x_worker_B'], e['driver_maxrss_B'], e['driver_fits_in_formula_headroom'], e['corrected_memory_max_B(1.3 x (driver + W x per_worker), page_up)'], e['raise_once_bytes'], kind))
    print('DECLARE: ' + out['declare'])
    return 0


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    raise SystemExit(main())

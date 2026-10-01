"""boxwalk.driver — prime-batched walk + farm produce() with gate battery.

walk():   replay the plan for a batch of K primes at once (windows stacked
          on a leading prime axis, the M-rail pattern: the recorded
          program is prime-independent, so one python pass serves K primes).
produce(): multiprocessing farm over prime batches -> CRT -> Wang rational
          reconstruction -> exact Fraction, gated:
  G1 crossover: walk vs dense oracle at a truncated u (value AND window).
  G2 mutation: corrupt one recorded coefficient -> value must change /
     replay must fail loudly.
  G3 two disjoint CRT sets reconstruct byte-identical Fractions.
Receipt manifest: spec/program/generator hashes + per-prime residues.
Logs: time-stamped per-step lines (rate-fittable, tools/eta.py style).
"""
import os, json, time, copy
from fractions import Fraction
import numpy as np

from . import core
from .relations import refill_batch, DegeneracyError

INT = np.int64


def _log_line(fh, msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
    if fh:
        fh.write(line + '\n')
        fh.flush()


def walk(spec, plan, primes, step_cb=None, log_fh=None):
    """Replay the plan for primes (list). Returns
    {'Z': [int mod p_i], 'u': final u, 'valid': per-axis valid extent}.
    step_cb(u, M, valid, pvec) is called after every step (selftest hook).
    """
    pvec = np.array([int(p) for p in primes], dtype=INT)
    S, n = plan['side'], plan['n']
    M = core.init_window_batch(S, n, pvec)
    valid = [S] * n
    u = {}
    nstep = 0
    for seg in plan['segments']:
        c = seg['class']
        poly = spec.polys[c]
        for step in range(seg['steps']):
            t0 = time.time()
            u = {**u, c: u.get(c, 0) + 1}
            Mnew, newvalid = core.raise_batch(M, poly, pvec, valid)
            if seg['mech'] == 'refill':
                refill_batch(Mnew, M, seg['program'], seg['tables'],
                             u, pvec)
                newvalid = [S] * n
            M, valid = Mnew, newvalid
            nstep += 1
            _log_line(log_fh,
                      f"boxwalk step={nstep} seg={c} i={step+1}/"
                      f"{seg['steps']} N={sum(u.values())} K={len(pvec)} "
                      f"dt={time.time()-t0:.3f}s")
            if step_cb:
                step_cb(u, M, valid, pvec)
    Z = M[(slice(None),) + (0,) * n] % pvec
    return {'Z': [int(z) for z in Z], 'u': dict(u), 'valid': list(valid)}


def _slab_only_plan(spec, u, n, min_side=3):
    """Ad-hoc all-slab plan for a small target u (no relations needed)."""
    order = sorted(u, key=lambda k: (u[k], k))
    S0 = max(min_side, 1 + sum(u[k] * spec.maxdeg(k) for k in order))
    return {'name': spec.name + '_trunc', 'spec_sha': spec.sha(),
            'order': order, 'side': S0, 'n': n, 'program_sha': 'slab-only',
            'segments': [{'class': k, 'steps': u[k], 'mech': 'slab',
                          'degs': spec.degs(k)} for k in order]}


def gate_crossover(spec, plan, primes=None, trunc_cap=2, log=print):
    """G1: dense-oracle crossover at truncated u (value + valid window),
    2 primes. Uses a slab-only walk (raise machinery + init window +
    mod-p arithmetic end-to-end)."""
    primes = primes or core.primes31(2)
    u_t = {k: min(v, trunc_cap) for k, v in spec.target_u.items()}
    pl = _slab_only_plan(spec, u_t, spec.n)
    res = walk(spec, pl, primes)
    ok = True
    for i, p in enumerate(primes):
        zo = core.Z_mod(spec, u_t, p)
        ok &= (res['Z'][i] % p == zo)
    # full-window compare, one prime, capture final window
    holder = {}

    def cb(u, M, valid, pvec):
        holder['M'], holder['valid'] = M, list(valid)
    res2 = walk(spec, pl, [primes[0]], step_cb=cb)
    v = min(holder['valid'])
    sub = holder['M'][0][tuple(slice(0, v) for _ in range(spec.n))]
    ref = core.moment_window(spec, u_t, v - 1, primes[0])
    ok &= bool(np.array_equal(sub % primes[0], ref % primes[0]))
    log(f"[gate G1 crossover] u_trunc={u_t} pass={bool(ok)}")
    return {'pass': bool(ok), 'u_trunc': u_t, 'primes': primes}


def gate_mutation(spec, plan, prime=None, log=print):
    """G2: corrupt one recorded coefficient (refill table if any, else a
    polynomial coefficient) -> Z must change or replay must fail loudly."""
    prime = prime or core.primes31(1)[0]
    detected = False
    refill_segs = [i for i, s in enumerate(plan['segments'])
                   if s['mech'] == 'refill']
    if refill_segs:
        # compare the FULL window at the end of the first refill segment
        # (a corrupted rim cell need not propagate to the corner Z)
        iseg = refill_segs[0]
        nstop = sum(s['steps'] for s in plan['segments'][:iseg + 1])
        seg = plan['segments'][iseg]

        def capture(pl_):
            snaps = []

            def cb(u, M, valid, pvec):
                if len(snaps) < nstop:
                    snaps.append(M.copy())
            walk(spec, pl_, [prime], step_cb=cb)
            return snaps[nstop - 1]
        M0 = capture(plan)
        plan2 = copy.deepcopy(plan)
        ti0 = seg['program'][0][0]
        t = plan2['segments'][iseg]['tables'][ti0]
        k0 = t.offsets[len(t.offsets) // 2]
        t.A[k0] = t.A[k0] + 12345
        try:
            M1 = capture(plan2)
            detected = not np.array_equal(M0, M1)
        except (DegeneracyError, AssertionError):
            detected = True
    else:
        z0 = walk(spec, plan, [prime])['Z'][0]
        spec2 = core.Spec(spec.to_dict())
        lab = plan['order'][0]
        k0 = next(iter(spec2.polys[lab]))
        spec2.polys[lab][k0] += 1
        z1 = walk(spec2, plan, [prime])['Z'][0]
        detected = (z1 != z0)
    log(f"[gate G2 mutation] corrupted-{'program' if refill_segs else 'poly'}"
        f" detected={detected}")
    return {'pass': bool(detected), 'mode':
            'program' if refill_segs else 'poly'}


def _worker(args):
    """Farm worker: one prime batch. Top-level for pickling."""
    spec_dict, plan, primes, logdir = args
    spec = core.Spec(spec_dict)
    fh = None
    if logdir:
        fh = open(os.path.join(
            logdir, f'walk_p{primes[0]}_K{len(primes)}.log'), 'a')
    try:
        res = walk(spec, plan, primes, log_fh=fh)
        return {'primes': primes, 'Z': res['Z'], 'u': res['u'], 'ok': True}
    except DegeneracyError as ex:
        return {'primes': primes, 'ok': False, 'error': str(ex)}
    finally:
        if fh:
            fh.close()


def produce(spec, plan=None, nprimes=8, batch=4, procs=None, outdir=None,
            gates=True, log=print, **plan_kw):
    """Exact Z(target_u) as a Fraction via nprimes-prime CRT + Wang ratrec,
    with the gate battery. Farm = multiprocessing pool over prime batches
    (launchd-compatible: pure python entry, per-step time-stamped logs).
    Returns manifest dict (also written to outdir/MANIFEST.json)."""
    from multiprocessing import Pool
    t00 = time.time()
    if plan is None:
        from . import planner
        plan = planner.plan(spec, log=log, **plan_kw)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    manifest = {'tool': 'boxwalk', 'name': spec.name,
                'spec_sha': spec.sha(), 'program_sha': plan['program_sha'],
                'target_u': spec.target_u, 'order': plan['order'],
                'side': plan['side'], 'nprimes': nprimes, 'batch': batch,
                'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                'gates': {}, 'ok': False}
    # ---- gates first (cheap, loud) --------------------------------------
    if gates:
        g1 = gate_crossover(spec, plan, log=log)
        g2 = gate_mutation(spec, plan, log=log)
        manifest['gates']['G1_crossover'] = g1
        manifest['gates']['G2_mutation'] = g2
        if not (g1['pass'] and g2['pass']):
            manifest['error'] = 'gate failure'
            _dump(manifest, outdir)
            return manifest
    # ---- farm ------------------------------------------------------------
    primes = core.primes31(nprimes)
    batches = [primes[i:i + batch] for i in range(0, nprimes, batch)]
    jobs = [(spec.to_dict(), plan, b, outdir) for b in batches]
    if procs is None:
        procs = min(len(batches), os.cpu_count() or 2)
    if procs > 1 and len(batches) > 1:
        with Pool(procs) as pool:
            results = pool.map(_worker, jobs)
    else:
        results = [_worker(j) for j in jobs]
    res_p, res_z = [], []
    for r in results:
        if not r['ok']:
            log(f"[produce] WARNING dropped batch {r['primes'][0]}..: "
                f"{r['error']}")
            continue
        res_p += r['primes']
        res_z += r['Z']
    manifest['n_primes_used'] = len(res_p)
    manifest['residues'] = {str(p): z for p, z in zip(res_p, res_z)}
    if len(res_p) < 4:
        manifest['error'] = 'too few surviving primes'
        _dump(manifest, outdir)
        return manifest
    # ---- G3: two disjoint CRT sets byte-identical ------------------------
    A = list(range(0, len(res_p), 2))
    B = list(range(1, len(res_p), 2))
    frs = []
    for idx in (A, B):
        r, M = core.crt_list([res_z[i] for i in idx],
                             [res_p[i] for i in idx])
        frs.append(core.ratrec(r, M))
    g3 = (frs[0] is not None and frs[1] is not None
          and repr(frs[0]).encode() == repr(frs[1]).encode())
    manifest['gates']['G3_disjoint_crt'] = {
        'pass': bool(g3),
        'half_A': str(frs[0]), 'half_B': str(frs[1])}
    log(f"[gate G3 disjoint-CRT byte-identical] pass={bool(g3)}")
    if not g3:
        manifest['error'] = ('disjoint CRT halves disagree or ratrec failed '
                             '— increase nprimes')
        _dump(manifest, outdir)
        return manifest
    # final value from ALL primes; must reduce to every residue
    r, M = core.crt_list(res_z, res_p)
    fr = core.ratrec(r, M)
    check = all(fr.numerator % p * pow(fr.denominator % p, p - 2, p) % p
                == z for p, z in zip(res_p, res_z))
    manifest['ok'] = bool(check)
    manifest['Z'] = {'num': str(fr.numerator), 'den': str(fr.denominator)}
    manifest['wall_s'] = round(time.time() - t00, 1)
    _dump(manifest, outdir)
    log(f"[produce] Z = {str(fr)[:80]}{'...' if len(str(fr)) > 80 else ''} "
        f"({manifest['wall_s']}s, {len(res_p)} primes)")
    manifest['Z_fraction'] = fr
    return manifest


def _dump(manifest, outdir):
    if outdir:
        m = {k: v for k, v in manifest.items() if k != 'Z_fraction'}
        json.dump(m, open(os.path.join(outdir, 'MANIFEST.json'), 'w'),
                  indent=1)

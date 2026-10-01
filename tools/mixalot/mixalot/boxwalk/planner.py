"""boxwalk.planner — class-agnostic segment planner.

Order: user-supplied list, or the 'ascending-tangency greedy' default:
walk classes in ascending target count (the largest class LAST, as the final
pure ray — the pilot's validated order rule generalized: early segments keep
the settled tangency set small, where relation modules are richest).

Per segment (class c, settled set C):
  - measure usable directions with a peel probe (pilot p3d/p3b machinery):
    segment 1 (C empty) uses free fields (validated full-cube refill);
    later segments use certified tangency-module generators (Singular syz,
    sampled-scan fallback).
  - mechanism = 'refill' iff the probe CLOSES the peel (all axes refill);
    else 'slab' (raise-only; the window is pre-grown and consumed at
    deg_e(P_c) layers/axis/step — pilot verdict: multi-class modules leak,
    slab the dead directions).

Window side S0 = base + sum over slab segments of steps*maxdeg(class).
Memory (per prime-batch job): K * 2 levels * S0^n * 8 bytes.
"""
import time
import numpy as np

from . import core, syz
from .relations import (RelationSystem, free_families, build_tables,
                        record_program)


def default_order(spec):
    """Ascending target count; ties broken by label. Largest class last."""
    return sorted(spec.target_u, key=lambda k: (spec.target_u[k], k))


def plan(spec, ref_prime=None, probe_side=None, min_side=6,
         syz_timeout=600, try_refill=True, log=print):
    """Build a walk plan: segment order, mechanisms, recorded programs.
    Programs are recorded at the production side S0, reference prime, and
    the u of each segment's first step (coefficients are affine in u; the
    replay asserts pivot non-vanishing at every later (u, prime))."""
    t00 = time.time()
    ref_prime = ref_prime or core.primes31(1)[0]
    order = (list(spec.order) if isinstance(spec.order, list)
             else default_order(spec))
    assert set(order) == set(spec.target_u), \
        f"order {order} != target classes {sorted(spec.target_u)}"
    n = spec.n
    system = RelationSystem(spec)

    # ---- pass 1: probe mechanisms at a small side --------------------
    pside = probe_side or max(5, min(7, min_side))
    seg_meta = []
    settled = []
    for i, c in enumerate(order):
        steps = spec.target_u[c]
        info = {'class': c, 'steps': steps, 'mech': 'slab',
                'degs': spec.degs(c), 'gen_meta': None, 'probe': None}
        u_probe = {k: spec.target_u[k] for k in settled}
        u_probe[c] = 1
        fams = None
        if try_refill and i == 0:
            fams = free_families(system, 2, [c])
            info['gen_meta'] = {'source': 'free_fields', 'n_certified':
                                len(fams)}
        elif try_refill:
            fams, gmeta = syz.certified_families(
                spec, system, settled, c, timeout=syz_timeout, log=log)
            info['gen_meta'] = gmeta
        if fams:
            tabs = build_tables(fams)
            known0 = np.zeros((pside,) * n, dtype=bool)
            core_sl = tuple(slice(0, pside - d) for d in info['degs'])
            known0[core_sl] = True
            stats, _, _ = record_program(tabs, u_probe, pside, ref_prime, c,
                                         down_ok=frozenset({c}),
                                         known0=known0)
            info['probe'] = stats
            if stats['closed']:
                info['mech'] = 'refill'
                info['_fams'] = fams
        log(f"[plan] seg{i} class={c} steps={steps} mech={info['mech']} "
            f"probe={info['probe'] and info['probe']['axis_face_coverage']}")
        seg_meta.append(info)
        settled.append(c)

    # ---- window size ---------------------------------------------------
    slab_budget = sum(s['steps'] * max(s['degs']) for s in seg_meta
                      if s['mech'] == 'slab')
    S0 = max(min_side, 1 + slab_budget)
    mem_gb = 2 * S0 ** n * 8 / 1e9      # per prime, 2 levels
    log(f"[plan] side S0={S0} (slab budget {slab_budget}); "
        f"memory {mem_gb:.3f} GB/prime (x K primes per batch job)")

    # ---- pass 2: record production programs at side S0 ------------------
    segments = []
    settled = []
    for i, s in enumerate(seg_meta):
        c = s['class']
        seg = {k: v for k, v in s.items() if k not in ('_fams',)}
        if s['mech'] == 'refill':
            u0 = {k: spec.target_u[k] for k in settled}
            u0[c] = 1
            tabs = build_tables(s['_fams'])
            known0 = np.zeros((S0,) * n, dtype=bool)
            core_sl = tuple(slice(0, S0 - d) for d in s['degs'])
            known0[core_sl] = True
            t0 = time.time()
            stats, program, usable = record_program(
                tabs, u0, S0, ref_prime, c, down_ok=frozenset({c}),
                known0=known0)
            if not stats['closed']:
                log(f"[plan] seg{i} {c}: refill did NOT close at S0={S0} "
                    f"({stats['n_unknown']} unknown) -> demote to slab")
                seg['mech'] = 'slab'
                seg['probe_S0'] = stats
            else:
                seg['program'] = program
                seg['tables'] = usable
                seg['probe_S0'] = stats
                log(f"[plan] seg{i} {c}: program recorded "
                    f"({stats['program_groups']} groups, "
                    f"{time.time()-t0:.1f}s)")
        segments.append(seg)
        settled.append(c)

    # slab budget may have grown if a refill was demoted: recheck
    slab_budget2 = sum(s['steps'] * max(s['degs']) for s in segments
                      if s['mech'] == 'slab')
    if 1 + slab_budget2 > S0:
        log("[plan] slab budget grew after demotion; replanning with "
            "try_refill unchanged, larger side")
        return plan(spec, ref_prime, probe_side, 1 + slab_budget2,
                    syz_timeout, try_refill, log)

    pl = {'name': spec.name, 'spec_sha': spec.sha(), 'order': order,
          'segments': segments, 'side': S0, 'n': n,
          'ref_prime': ref_prime, 'mem_gb_per_prime': mem_gb,
          'wall_s': round(time.time() - t00, 1)}
    pl['program_sha'] = program_sha(pl)
    return pl


def program_sha(pl):
    """Receipt-style witness hash of the recorded pivot programs +
    generator tables (tools/trust/receipt design): the program + certified
    generators ARE a per-(u-path, prime) certificate that every refilled
    value is a linear consequence of certified boundary-safe relations."""
    import hashlib
    h = hashlib.sha256()
    h.update(pl['spec_sha'].encode())
    for seg in pl['segments']:
        h.update(seg['class'].encode())
        h.update(seg['mech'].encode())
        for ti, dstar, pos in seg.get('program', []):
            h.update(bytes(f"{ti}|{dstar}", 'ascii'))
            h.update(np.ascontiguousarray(pos).tobytes())
        for t in seg.get('tables', []):
            h.update(repr(sorted(t.A.items())).encode())
            h.update(repr(sorted((k, tuple(v)) for k, v in
                                 t.B.items())).encode())
            h.update(repr(sorted((k, sorted(v.items())) for k, v in
                                 t.C.items())).encode())
            h.update(repr(sorted((k, sorted(v.items())) for k, v in
                                 t.down.items())).encode())
    return h.hexdigest()[:16]

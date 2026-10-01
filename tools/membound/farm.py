r"""membound.farm — process farm for omega-core quadratures.

Every unique (core, eps-order, sign-domain class) quadrature of
membound.core.compute_core is an independent task; this module farms them
over a process pool and re-assembles the eps-layer moments. It is the
parallel path behind `python3 -m membound --spec ... --jobs N` and behind
the c_M validation gate (gate_cM.py). Each worker rebuilds the core from its
plain-JSON form (OmegaCoreSpec.to_dict / from_dict), so any user-supplied
core farms the same way the registry cores do; nothing is hard-wired.
"""
import concurrent.futures as cf
import time

import mpmath as mp

from .core import OmegaCoreSpec, KCache, make_integrand, domain_classes


def parse_panels(panels):
    """'0,0.5,2,8,inf' (or a list of such tokens) -> mpf breakpoints."""
    toks = panels.split(",") if isinstance(panels, str) else list(panels)
    toks = [str(t).strip() for t in toks]
    brk = [mp.inf if t == "inf" else mp.mpf(t) for t in toks]
    if len(brk) < 2:
        raise ValueError("need at least two panel breakpoints, e.g. 0,inf")
    return toks, brk


def _quad_task(task):
    """Worker: one (core, order, domain-class) quad with its own dps + KCache."""
    mp.mp.dps = task["dps"]
    spec = OmegaCoreSpec.from_dict(task["spec"])
    kc = KCache(task["dps"])
    _, brk = parse_panels(task["brk"])
    rep = tuple(task["dom"])
    f = make_integrand(spec, rep, task["order"], kc)
    t0 = time.time()
    if spec.n_freq == 3:
        third = [mp.mpf(0), mp.mpf(1)] if rep[0] == 'gt' else brk
        val = mp.quad(f, brk, brk, third, maxdegree=task["maxdegree"])
    else:
        val = mp.quad(f, brk, brk, maxdegree=task["maxdegree"])
    return (task["key"], task["order"], task["mult"], repr_num(val),
            time.time() - t0)


def repr_num(v):
    """Loss-free string form of an mpf/mpc at the working dps (worker -> parent)."""
    if isinstance(v, mp.mpc):
        return ("c", mp.nstr(v.real, mp.mp.dps + 5), mp.nstr(v.imag, mp.mp.dps + 5))
    return ("r", mp.nstr(v, mp.mp.dps + 5))


def parse_num(t):
    if t[0] == "c":
        return mp.mpc(mp.mpf(t[1]), mp.mpf(t[2]))
    return mp.mpf(t[1])


def compute_cores_parallel(specs, dps, panels, maxdegree, jobs, log=print,
                           orders=(0, 1, 2)):
    """Farm all unique (core, order, class) quads of several cores at once.

    specs : dict key -> OmegaCoreSpec
    panels: '0,0.5,2,8,inf' or token list (magnitude-axis breakpoints)
    Returns dict key -> {order: J^(order)} INCLUDING the 1/(2pi)^n measure,
    i.e. the same object compute_core returns per core.
    """
    mp.mp.dps = dps
    toks, _ = parse_panels(panels)
    kc = KCache(dps)
    tasks = []
    for key, spec in specs.items():
        for order in orders:
            for rep, mult, _ in domain_classes(spec, order, kc):
                tasks.append({"key": key, "spec": spec.to_dict(), "order": order,
                              "dom": list(rep), "mult": mult, "dps": dps,
                              "brk": toks, "maxdegree": maxdegree})
    log(f"  farming {len(tasks)} quads across {jobs} processes")
    acc = {key: {o: mp.mpf(0) for o in orders} for key in specs}
    if jobs <= 1:
        results = map(_quad_task, tasks)
    else:
        ex = cf.ProcessPoolExecutor(max_workers=jobs)
        results = ex.map(_quad_task, tasks)
    try:
        for key, order, mult, tval, dt in results:
            val = parse_num(tval)
            acc[key][order] += mult * val
            log(f"    [{specs[key].name} k={order}] x{mult} "
                f"quad={mp.nstr(val, min(20, dps))} ({dt:.1f}s)")
    finally:
        if jobs > 1:
            ex.shutdown()
    out = {}
    for key, spec in specs.items():
        inv = 1 / (2 * mp.pi) ** spec.n_freq
        out[key] = {o: inv * v for o, v in acc[key].items()}
    return out

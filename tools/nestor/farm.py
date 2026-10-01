#!/usr/bin/env python3
r"""nestor.farm -- node-farmed explicit tanh-sinh grid (the oracle_par engine).

THE FARM ENGINE (verbatim ported math; NOT rewritten)
=====================================================
  (measured: C_2222(21i) dps=80 md=8 = 2573 nodes,
        360 s wall at nproc=48 on a heavily loaded host, ~100x vs serial; smoke
        C_22 = E4/pi^4 to 44.3 d; farmed C_2222(2i) = serial mp.quad to all
        43 digits).
  The originating oracle_par.py (not part of this repo) has its worker
  layer welded to one oracle's H(u) dispatch, so this module PORTS the
  engine and generalizes ONLY the plumbing (integrand by importable spec,
  interval [a,b] instead of the hard-wired [0,1/2]); every
  grid/weight/nesting formula is the source's.

THE ENGINE (source math, general [a,b])
=======================================
  Outer map (oracle_par.py header, u = (1+x)/4 generalized to affine):
    x = tanh((pi/2) sinh t), t = k*h, h = 2^-md, s = (pi/2) sinh t;
    cancellation-free fraction  (1+x)/2 = 1/(1+e^{-2s})
    (source: u = 1/(2(1+e^{-2s})) IS this with [a,b]=[0,1/2]);
    u_k = a + (b-a)/(1+e^{-2s_k}),
    W_k = (b-a)/2 * h * (pi/2) cosh(t_k)/cosh^2(s_k),
    int_a^b f du = sum_k W_k f(u_k).
    (source normalization check: (b-a)/2 = 1/4 on [0,1/2] reproduces
     "int_0^{1/2} H du = (1/4) sum_k w_k H(u_k)" exactly.)
  Node cutoff (_nodes, verbatim): symmetric -K..K, stop when
    h*(pi/2)cosh(t)/cosh^2(s) * f_bound < 10^-(dps+15).
  Nesting (verbatim): tanh-sinh levels nest -- level md-1 = the even nodes
    of level md -- so one farmed grid yields the md-1 <-> md
    self-consistency check for free:
      I_hi = sum_all W_k f_k     (level md)
      I_lo = 2 * sum_{k even} W_k f_k   (level md-1: h doubles => weights x2)
    selfcons digits = -log10 |I_hi - I_lo| / |I_hi|.
  Farming (verbatim pattern): chunks ks[i::nch], nch = max(nproc*6, 1),
  Pool(initializer=...) + imap_unordered, values transported as
  mp.nstr(dps+15) strings, per-chunk live rate + ETA log (timed-pilot rate-fit),
  per-config JSON checkpoint with resume = skip (receipts.checkpoint_book).

SWEEP-PROOF POOL DISCIPLINE (pool-parent liveness; no orphan pools)
===================================================================
  * The Pool lives in a with-block AND a try/finally terminate/join.
  * Every worker records the parent pid at init and aborts its chunk loop
    if it finds itself reparented (os.getppid() changed) -- a SIGKILLed
    parent cannot leave a silently grinding orphan pool.
  * The parent pid + worker pids go into the receipt (the kill handle).

"""
from __future__ import annotations

import importlib
import multiprocessing as mproc
import os
import sys
import time

import mpmath as mp
from mpmath import mpf

from . import receipts as _receipts
from .ladder import MissingRefusal

_G = {}


def _load_spec(spec):
    """Resolve 'module:function' (or a callable) to a callable.
    Workers re-import by spec -- the source pattern (_winit imports
    its oracle module by name); a bare callable is only usable serially."""
    if callable(spec):
        return spec
    if not isinstance(spec, str) or ":" not in spec:
        raise MissingRefusal(
            "farm integrand spec must be 'module:function' (importable in "
            "workers) or a callable (serial only); got %r" % (spec,))
    modname, fname = spec.split(":", 1)
    try:
        mod = importlib.import_module(modname)
    except ImportError as e:
        # unimportable spec = missing dependency (exit 4), NOT an
        # evaluation crash (a raw escape here would be renamed
        # EvalRefusal/exit 3 by the front door — the wrong class of the
        # contract)
        raise MissingRefusal(
            "farm integrand spec %r: module %s not importable in this "
            "process (%s: %s)" % (spec, modname, type(e).__name__, e))
    f = getattr(mod, fname, None)
    if f is None:
        raise MissingRefusal("farm integrand spec %r: %s has no attribute %s"
                             % (spec, modname, fname))
    return f


def _winit(spec, a_str, b_str, dps, md, path_entry, ppid, fkw, wants_dists,
           in_worker=False):
    """Worker initializer (ported pattern: oracle_par._winit)."""
    import mpmath as _mp
    if in_worker:
        # pool workers must be PLAINLY killable: the pool parent installs a
        # SIGTERM->AbortRefusal handler (oracle.py abort discipline) and
        # forked workers inherit it -- inside a worker that handler turns
        # pool.terminate()'s SIGTERM into a swallowed task exception and
        # pool.join() deadlocks here (measured) — use the reset below.
        # Reset to the default disposition IN WORKERS ONLY (serial mode runs
        # _winit in the parent process, whose abort handler must survive);
        # the getppid orphan guard and pool.terminate() stay the kill
        # handles.
        import signal as _sig
        try:
            _sig.signal(_sig.SIGTERM, _sig.SIG_DFL)
        except (ValueError, OSError):          # non-main thread: leave it
            pass
    _mp.mp.dps = dps + 20                      # source: dps + 20
    if path_entry and path_entry not in sys.path:
        sys.path.insert(0, path_entry)
    _G.update(mp=_mp, f=_load_spec(spec), a=_mp.mpf(a_str), b=_mp.mpf(b_str),
              dps=dps, h=_mp.mpf(2) ** (-md), ppid=ppid, fkw=fkw or {},
              wants_dists=bool(wants_dists))


def _node_uw(_mp, k, h, a, b):
    """(u_k, d_a, d_b, W_k) by the cancellation-free source formulas.
    d_a = u-a and d_b = b-u are BOTH computed to full RELATIVE accuracy
    (the complement-pair form of the originating evaluator's w_ts_nodes, which returns
    (xm, xm_c, w) for exactly this reason; oracle_par's u = 1/(2(1+e^-2s))
    is the d_a formula on [0,1/2])."""
    t = k * h
    s = (_mp.pi / 2) * _mp.sinh(t)
    ex = _mp.exp(-2 * s)
    d_a = (b - a) / (1 + ex)                   # u - a, relative-accurate
    d_b = (b - a) * ex / (1 + ex)              # b - u, relative-accurate
    w = (b - a) / 2 * h * (_mp.pi / 2) * _mp.cosh(t) / _mp.cosh(s) ** 2
    return a + d_a, d_a, d_b, w


def _wchunk(ks):
    """Evaluate W_k f(u_k) over one chunk (ported pattern: _wchunk);
    returns string-transported values + chunk wall.  Orphan guard: abort
    if the pool parent is gone (reparented)."""
    _mp = _G['mp']
    f = _G['f']
    a, b, dps, h = _G['a'], _G['b'], _G['dps'], _G['h']
    fkw = _G['fkw']
    wd = _G['wants_dists']
    out = []
    t0 = time.time()
    for k in ks:
        if os.getppid() != _G['ppid']:         # sweep-proof: no orphan pools
            raise SystemExit(3)
        u, d_a, d_b, w = _node_uw(_mp, k, h, a, b)
        v = f(u, d_a, d_b, **fkw) if wd else f(u, **fkw)
        out.append((k, _mp.nstr(w * v, dps + 15, strip_zeros=False)))
    return out, time.time() - t0


def ts_nodes(a, b, dps, md, f_bound, exp_a=0, exp_b=0):
    """Node list -K_neg..K_pos (ported from oracle_par._nodes, interval
    generalized): each side stops when weight x |f|-bound < 10^-(dps+15)
    -- the source cutoff law.  For a spec-declared ALGEBRAIC endpoint
    singularity the bound at a node is f_bound * d_a^exp_a * d_b^exp_b
    (d_a, d_b cancellation-free; exp = 0 reproduces the source criterion
    exactly and gives a symmetric -K..K list)."""
    a, b = mpf(a), mpf(b)
    h = mpf(2) ** (-md)
    Hb = mpf(f_bound)
    tol = mpf(10) ** (-(dps + 15))
    ea, eb = mpf(exp_a), mpf(exp_b)
    Ks = []
    for sgn in (1, -1):                        # +K side: d_b small; -K: d_a
        K = 0
        while True:
            K += 1
            _u, d_a, d_b, w = _node_uw(mp, sgn * K, h, a, b)
            # source cutoff law tests the RAW map weight (no interval
            # factor): h*(pi/2)cosh t/cosh^2 s * bound < tol -- byte-parity
            # of the node SET with oracle_par._nodes requires it verbatim
            # (H1 parity check)
            w_raw = w / ((b - a) / 2)
            bound = Hb
            if ea:
                bound *= d_a ** ea
            if eb:
                bound *= d_b ** eb
            if w_raw * bound < tol:
                break
        Ks.append(K)
    return list(range(-Ks[1], Ks[0] + 1))


def farm_integrate(spec, a, b, dps, md, nproc, f_bound=1, f_kwargs=None,
                   exp_a=0, exp_b=0, wants_dists=False,
                   path_entry=None, log=print, tag="farm"):
    """Node-farmed tanh-sinh integral of spec over [a,b] at level md
    (ported pattern: oracle_par.eval_point).

    spec       'module:function' importable in workers (f(u, **f_kwargs) ->
               mpf), or a callable (then nproc is forced to 1: bare
               callables do not cross Pool pickling -- stated, not hidden).
    f_bound    bound on |f| for the node cutoff (source: Hb).
    exp_a/exp_b declared algebraic endpoint exponents (spec-verified by
               singcheck); enter the node-cutoff bound as
               f_bound * d_a^exp_a * d_b^exp_b.
    wants_dists integrand signature f(u, d_a, d_b, **kw): receives the
               cancellation-free endpoint distances (needed for full
               relative accuracy of a singular factor near an endpoint).
    path_entry extra sys.path entry for worker imports (computed by the
               caller at runtime; nestor itself hard-codes NO paths).
    Returns a record dict (the oracle_par checkpoint schema, extended):
      value/value_lo (nstr dps+10), selfcons_digits, nodes, wall_s,
      a/b/dps/md/nproc, pool pids (kill handle).
    """
    old = mp.mp.dps
    mp.mp.dps = dps + 20                       # source: dps + 20
    try:
        a_str, b_str = str(a), str(b)
        ks = ts_nodes(a, b, dps, md, f_bound, exp_a=exp_a, exp_b=exp_b)
        if callable(spec) and nproc > 1:
            log("[%s] bare-callable spec: forcing nproc=1 (serial); pass "
                "'module:function' to farm" % tag)
            nproc = 1
        nch = max(nproc * 6, 1)                # source: nproc*6 chunking
        chunks = [ks[i::nch] for i in range(nch)]
        chunks = [c for c in chunks if c]
        log("[%s] int_[%s,%s] dps=%d md=%d nodes=%d chunks=%d nproc=%d"
            % (tag, a_str, b_str, dps, md, len(ks), len(chunks), nproc))
        t0 = time.time()
        vals = {}
        done = 0
        pids = []
        if nproc <= 1:
            # serial: no pool, so the orphan guard must be a no-op --
            # _wchunk compares os.getppid(), which in-process is OUR parent
            _winit(spec, a_str, b_str, dps, md, path_entry, os.getppid(),
                   f_kwargs, wants_dists)
            for c in chunks:
                res, _secs = _wchunk(c)
                for k, vstr in res:
                    vals[k] = mpf(vstr)
                done += len(res)
        else:
            pool = mproc.Pool(nproc, initializer=_winit,
                              initargs=(spec, a_str, b_str, dps, md,
                                        path_entry, os.getpid(), f_kwargs,
                                        wants_dists, True))
            try:
                with pool:
                    pids = [p.pid for p in pool._pool]
                    for res, _secs in pool.imap_unordered(_wchunk, chunks):
                        for k, vstr in res:
                            vals[k] = mpf(vstr)
                        done += len(res)
                        el = time.time() - t0
                        rate = done / el
                        eta = (len(ks) - done) / rate if rate > 0 \
                            else float('inf')
                        log("  %d/%d nodes  %7.1fs  %6.2f nodes/s  eta %6.0fs"
                            % (done, len(ks), el, rate, eta))
            finally:
                pool.terminate()               # no orphan pools
                pool.join()
        # nesting: level md-1 = even nodes; general-weight form (see header)
        S = mp.fsum(vals.values())
        S_even = mp.fsum(v for k, v in vals.items() if k % 2 == 0)
        I_hi = S                               # level md
        I_lo = 2 * S_even                      # level md-1 (h doubles)
        sc = (float(-mp.log10(abs(I_hi - I_lo) / abs(I_hi)))
              if I_hi != I_lo else 999.0)      # source selfcons measure
        secs = time.time() - t0
        log("  => %s  selfcons(md%d<->%d)=%.1fd  wall=%.0fs"
            % (mp.nstr(I_hi, min(dps, 40)), md - 1, md, sc, secs))
        return dict(a=a_str, b=b_str, dps=dps, md=md, nproc=nproc,
                    spec=(spec if isinstance(spec, str)
                          else getattr(spec, "__name__", "callable")),
                    value=mp.nstr(I_hi, dps + 10, strip_zeros=False),
                    value_lo=mp.nstr(I_lo, dps + 10, strip_zeros=False),
                    selfcons_digits=round(sc, 2), nodes=len(ks),
                    wall_s=round(secs, 1),
                    pool=dict(parent_pid=os.getpid(), worker_pids=pids))
    finally:
        mp.mp.dps = old


def farm_campaign(spec, points, out_path, log=print, **kw):
    """Checkpointed battery (ported pattern: oracle_par.main): one record
    per config, resume = skip, atomic tmp+replace writes via receipts."""
    book = _receipts.checkpoint_load(out_path)
    for p in points:
        key = "%s_%s_%s_%s_%s" % (spec, p["a"], p["b"], p["dps"], p["md"])
        if key in book:
            log("[skip] %s already done" % key)
            continue
        rec = farm_integrate(spec, p["a"], p["b"], p["dps"], p["md"],
                             p.get("nproc", 1), log=log, **kw)
        book[key] = rec
        _receipts.checkpoint_save(out_path, book)
        log("[ckpt] %s -> %s" % (key, out_path))
    return book

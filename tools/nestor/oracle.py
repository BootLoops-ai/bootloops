#!/usr/bin/env python3
r"""nestor.oracle -- the ONE front door: nestor.integrate / nestor.oracle.

NESTOR (NESted Tanh-sinh ORacle) evaluates
definite integrals to controlled high precision and reports ONLY the digits
successive refinements agree on.  This module wires the four legs of the
consolidated layer behind a single call:

  1. load-time singularity verification of the integrand spec  (singcheck)
  2. per-level precision ladder, refine-until-agreement, fail-closed (ladder)
  3. the node farm for expensive (typically nested) integrands       (farm)
  4. JSON receipts: per-level walls, agreement digits, refusals  (receipts)

ENTRY POINTS
============
  nestor.integrate(f, a, b, dps=30, ...) -> mpf
      The value at >= dps certified-by-agreement digits, or a NAMED
      refusal (NestorRefusal subclass; .exit_code follows the evaluators'
      0/2/3/4 contract).  A value never escapes a failed loop.
  nestor.oracle(f, a, b, dps=30, ...) -> receipt dict
      Same run, full receipt returned (value_str, rungs, singcheck record,
      refusal record on failure -- the refusal is ALSO raised unless
      raise_on_refusal=False).

ROUTES
======
  serial (default): cert_quad -- nested tanh-sinh depths escalated until
      the double-refinement agreement beats 10^-(dps+guard) (the
      originating certified-quadrature design, ported in ladder.py).
  farmed (nproc>1 or md set): the oracle_par.py explicit-grid engine
      (farm.py) -- one grid at level md farmed over a Pool, md-1<->md
      self-consistency free from grid nesting, PLUS a ladder over md rungs
      when the single grid's self-consistency misses the bar.

NESTING
=======
  An integrand may itself call nestor.integrate (the production nesting
  pattern: certification evaluators whose integrands integrate).
  Inner requests should use inner_tol_exp(D)/wp_for(D) from ladder.py so
  working precision scales per level -- the "wp-scaled ladder" of the
  NESTOR registration.

PROVENANCE: see the src-sha headers of ladder.py / farm.py; receipts stamp
the same shas into every run record.

"""
from __future__ import annotations

import signal as _signal
import time

import mpmath as mp

from . import farm as _farm
from . import receipts as _receipts
from . import singcheck as _singcheck
from .ladder import (AbortRefusal, CertFail, DEFAULT_GUARD, EvalRefusal,
                     NestorRefusal, cert_quad, digits, ladder, prec_for,
                     wp_for)


def oracle(f, a, b, dps=30, spec=None, guard=DEFAULT_GUARD, md=None,
           nproc=1, f_bound=None, f_kwargs=None, wants_dists=False,
           path_entry=None, receipt_path=None, tag="nestor", log=print,
           raise_on_refusal=True):
    """Full-receipt evaluation (see module docstring).  Returns the receipt
    dict; on refusal records it, writes it (if receipt_path), then re-raises
    unless raise_on_refusal=False.

    ABORT DISCIPLINE: while a run is active a
    SIGTERM is converted to AbortRefusal (recorded in the receipt, pool
    torn down by farm's finally, ALWAYS re-raised -- raise_on_refusal does
    not swallow aborts); the previous handler is restored on exit.  ANY
    other exception (an integrand crash like ZeroDivisionError at an
    interior pole node, KeyboardInterrupt, ...) is also recorded in the
    receipt before propagating -- a written receipt can never be silent
    about why there is no value.  Integrand crashes re-raise as the NAMED
    EvalRefusal with the original exception chained."""
    def _on_term(signum, frame):
        raise AbortRefusal(
            "nestor[%s]: SIGTERM received mid-run -- aborted fail-closed; "
            "partial rungs stay partial (no value)" % tag)
    # spec-declared endpoint singularity => the EXPLICIT-grid route even
    # serially: mpmath's TanhSinh rule saturates ~1e-32 on singular
    # endpoints (node cancellation) while the cancellation-free grid here
    # does not.
    singular = bool(spec) and any(
        spec.get(sd, {}).get("type", "none") != "none" for sd in ("a", "b"))
    farmed = (nproc and nproc > 1) or md is not None or singular \
        or wants_dists
    rec = _receipts.new_receipt(
        "farm" if farmed else "integrate",
        config=dict(a=str(a), b=str(b), dps=dps, guard=guard, md=md,
                    nproc=nproc, tag=tag,
                    spec_declared=spec is not None))
    t0 = time.time()
    installed, old_term = False, None
    try:
        try:
            import multiprocessing as _mproc
            # never install inside a daemon (pool-worker) process: workers
            # must stay plainly killable (see farm._winit; a nested nestor
            # call in a worker would otherwise reopen the terminate-deadlock
            # window for its duration)
            if not _mproc.current_process().daemon:
                old_term = _signal.signal(_signal.SIGTERM, _on_term)
                installed = True
        except ValueError:                     # not the main thread
            pass
        # ---- leg 1: load-time singularity verification (fail-closed) ----
        if isinstance(f, str) and path_entry:
            import sys
            if path_entry not in sys.path:
                sys.path.insert(0, path_entry)
        fx = _farm._load_spec(f) if isinstance(f, str) else f
        if wants_dists:
            fs = lambda x: fx(x, x - mp.mpf(a), mp.mpf(b) - x)
        else:
            fs = fx
        rec["singcheck"] = _singcheck.verify(fs, a, b, spec,
                                             wp=min(dps, 30))
        # ---- legs 2+3: ladder / farm ----
        if not farmed:
            wp = wp_for(dps, guard)
            with mp.workprec(prec_for(wp)):
                val, agr, depth = cert_quad(fx, a, b, dps + guard, wp,
                                            tag=tag)
                val = mp.re(val) if mp.im(val) == 0 else val
                _receipts.add_rung(rec, level=depth, dps=wp,
                                   wall_s=round(time.time() - t0, 2),
                                   value_str=mp.nstr(val, dps + 10,
                                                     strip_zeros=False),
                                   agreement_digits=round(
                                       float(-mp.log10(agr)) if agr > 0
                                       else float(wp), 2))
                if isinstance(val, mp.mpc):
                    # complex value: the receipt carries re/im as SEPARATE
                    # full-precision strings -- nstr's "(a + bj)" form is
                    # not re-parseable at precision (the originating
                    # evaluators keep complex through the folds and
                    # re-reduce at the end)
                    rec["value"] = mp.nstr(mp.re(val), dps + 10,
                                           strip_zeros=False)
                    rec["value_im"] = mp.nstr(mp.im(val), dps + 10,
                                              strip_zeros=False)
                else:
                    rec["value"] = mp.nstr(val, dps + 10, strip_zeros=False)
                rec["agreed_digits"] = rec["rungs"][-1]["agreement_digits"]
        else:
            md0 = md if md is not None else 6
            bound = f_bound if f_bound is not None else 1
            ea = eb = 0
            if singular:
                da = spec.get("a", {})
                db = spec.get("b", {})
                ea = float(da.get("exponent") or 0) \
                    if da.get("type") == "algebraic" else 0
                eb = float(db.get("exponent") or 0) \
                    if db.get("type") == "algebraic" else 0

            def eval_md(m):
                r = _farm.farm_integrate(f, a, b, dps, m, nproc,
                                         f_bound=bound, f_kwargs=f_kwargs,
                                         exp_a=ea, exp_b=eb,
                                         wants_dists=wants_dists,
                                         path_entry=path_entry, log=log,
                                         tag="%s.md%d" % (tag, m))
                _receipts.add_rung(rec, level=m, dps=dps,
                                   wall_s=r["wall_s"],
                                   value_str=r["value"],
                                   agreement_digits=r["selfcons_digits"],
                                   nodes=r["nodes"], pool=r["pool"])
                return mp.mpf(r["value"])

            # every value-string parse must happen at WORKING precision --
            # ambient-dps parsing truncates to ~15 dps
            with mp.workprec(prec_for(dps + 20)):
                # one grid may already certify: md-1<->md selfcons free
                v0 = eval_md(md0)
                if rec["rungs"][-1]["agreement_digits"] >= dps:
                    val, agr = v0, rec["rungs"][-1]["agreement_digits"]
                else:
                    # ladder over md rungs (refine-until-agreement,
                    # fail-closed)
                    val, agr, _ = ladder(eval_md,
                                         [md0 + 1, md0 + 2, md0 + 3],
                                         dps, tag=tag + ".mdladder",
                                         wp=dps + guard)
                rec["value"] = mp.nstr(val, dps + 10, strip_zeros=False)
            rec["agreed_digits"] = round(float(agr), 2)
        rec["wall_s"] = round(time.time() - t0, 2)
        return rec
    except NestorRefusal as e:
        _receipts.add_refusal(rec, e)
        rec["wall_s"] = round(time.time() - t0, 2)
        if raise_on_refusal or isinstance(e, AbortRefusal):
            raise                              # aborts are NEVER swallowed
        return rec
    except Exception as e:
        # integrand/machinery crash: record it, then re-raise NAMED with
        # the original chained (value-or-named-refusal entry contract)
        _receipts.add_refusal(rec, e)
        rec["refusal"]["named"] = False
        rec["wall_s"] = round(time.time() - t0, 2)
        raise EvalRefusal(
            "nestor[%s]: evaluation crashed (%s: %s) -- fail-closed, no "
            "value; original exception chained" %
            (tag, type(e).__name__, e)) from e
    except BaseException as e:
        # KeyboardInterrupt / SystemExit: record the abort, propagate as-is
        _receipts.add_refusal(rec, e)
        rec["refusal"]["named"] = False
        rec["wall_s"] = round(time.time() - t0, 2)
        raise
    finally:
        if installed:
            _signal.signal(_signal.SIGTERM, old_term)
        if receipt_path:
            _receipts.write(rec, receipt_path)


def integrate(f, a, b, dps=30, **kw):
    """The value-only front door: mpf at >= dps agreed digits, or a NAMED
    refusal (fail-closed -- no value escapes a failed loop)."""
    kw.pop("raise_on_refusal", None)
    rec = oracle(f, a, b, dps=dps, raise_on_refusal=True, **kw)
    # parse at working precision: an ambient-dps parse truncates to ~15 dps
    with mp.workprec(prec_for(dps + 20)):
        v = mp.mpf(rec["value"])
        vi = rec.get("value_im")
        return mp.mpc(v, mp.mpf(vi)) if vi is not None else v
